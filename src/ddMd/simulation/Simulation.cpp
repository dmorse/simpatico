#ifndef DDMD_SIMULATION_CPP
#define DDMD_SIMULATION_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Simulation.h"
#include <ddMd/storage/AtomIterator.h>
#include <ddMd/storage/GroupStorage.tpp>
#include <ddMd/integrators/Integrator.h>
#include <ddMd/integrators/IntegratorFactory.h>
#include <ddMd/configIos/ConfigIo.h>
#include <ddMd/configIos/ConfigIoFactory.h>
#include <ddMd/configIos/DdMdConfigIo.h>
#include <ddMd/configIos/SerializeConfigIo.h>
#include <ddMd/diagnostics/DiagnosticManager.h>

#ifndef DDMD_NOPAIR
#include <ddMd/potentials/pair/PairPotential.h>
#include <ddMd/potentials/pair/PairPotentialImpl.h>
#include <ddMd/potentials/pair/PairFactory.h>
#endif
#include <ddMd/potentials/bond/BondPotential.h>
#include <ddMd/potentials/bond/BondPotentialImpl.h>
#include <ddMd/potentials/bond/BondFactory.h>
#ifdef INTER_ANGLE
#include <ddMd/potentials/angle/AnglePotential.h>
#include <ddMd/potentials/angle/AngleFactory.h>
#endif
#ifdef INTER_DIHEDRAL
#include <ddMd/potentials/dihedral/DihedralPotential.h>
#include <ddMd/potentials/dihedral/DihedralFactory.h>
#endif
#ifdef INTER_EXTERNAL
#include <ddMd/potentials/external/ExternalPotential.h>
#include <ddMd/potentials/external/ExternalPotentialImpl.h>
#include <ddMd/potentials/external/ExternalFactory.h>
#endif

// namespace Util
#include <util/misc/FileMaster.h>
#include <util/ensembles/EnergyEnsemble.h>
#include <util/ensembles/BoundaryEnsemble.h>
#include <util/space/Vector.h>
#include <util/space/IntVector.h>
#include <util/space/Tensor.h>
#include <util/param/Factory.h>
#include <util/misc/Log.h>
#include <util/mpi/MpiSendRecv.h>
#include <util/misc/Timer.h>

#include <util/format/Int.h>
#include <util/format/Dbl.h>
#include <util/format/Str.h>
#include <util/misc/ioUtil.h>
#include <util/misc/initStatic.h>

// std headers
#include <fstream>
#include <unistd.h>
#include <stdlib.h>

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   #ifdef UTIL_MPI
   Simulation::Simulation(MPI::Intracomm& communicator)
   #else
   Simulation::Simulation()
   #endif
    : atomStorage_(),
      bondStorage_(),
      #ifdef INTER_ANGLE
      angleStorage_(),
      #endif
      #ifdef INTER_DIHEDRAL
      dihedralStorage_(),
      #endif
      boundary_(),
      atomTypes_(),
      domain_(),
      buffer_(),
      exchanger_(),
      random_(),
      maxBoundary_(),
      kineticEnergy_(0.0),
      pairPotentialPtr_(0),
      bondPotentialPtr_(0),
      #ifdef INTER_ANGLE
      anglePotentialPtr_(0),
      #endif
      #ifdef INTER_DIHEDRAL
      dihedralPotentialPtr_(0),
      #endif
      #ifdef INTER_EXTERNAL
      externalPotentialPtr_(0),
      #endif
      integratorPtr_(0),
      energyEnsemblePtr_(0),
      boundaryEnsemblePtr_(0),
      fileMasterPtr_(0),
      configIoPtr_(0),
      serializeConfigIoPtr_(0),
      diagnosticManagerPtr_(0),
      #ifndef DDMD_NOPAIR
      pairFactoryPtr_(0),
      #endif
      bondFactoryPtr_(0),
      #ifdef INTER_ANGLE
      angleFactoryPtr_(0),
      #endif
      #ifdef INTER_DIHEDRAL
      dihedralFactoryPtr_(0),
      #endif
      #ifdef INTER_EXTERNAL
      externalFactoryPtr_(0),
      #endif
      integratorFactoryPtr_(0),
      configIoFactoryPtr_(0),
      #ifndef DDMD_NOPAIR
      pairStyle_(),
      #endif
      bondStyle_(),
      #ifdef INTER_ANGLE
      angleStyle_(),
      #endif
      #ifdef INTER_DIHEDRAL
      dihedralStyle_(),
      #endif
      #ifdef INTER_EXTERNAL
      externalStyle_(),
      #endif
      nAtomType_(0),
      nBondType_(0),
      #ifdef INTER_ANGLE
      nAngleType_(0),
      #endif
      #ifdef INTER_DIHEDRAL
      nDihedralType_(0),
      #endif
      #ifdef INTER_EXTERNAL
      hasExternal_(false),
      #endif
      maskedPairPolicy_(MaskBonded),
      reverseUpdateFlag_(false),
      #ifdef UTIL_MPI
      communicator_(communicator),
      #endif
      isInitialized_(false),
      isRestarting_(false)
   {
      Util::initStatic();
      setClassName("Simulation");

      #ifdef UTIL_MPI
      if (!MPI::Is_initialized()) {
         UTIL_THROW("MPI is not initialized");
      }
      Vector::commitMpiType();
      IntVector::commitMpiType();
      AtomType::initStatic();

      setIoCommunicator(communicator);
      domain_.setGridCommunicator(communicator);
      #else
      domain_.setRank(0);
      #endif

      // Set connections between member objects
      domain_.setBoundary(boundary_);
      #if 0
      exchanger_.associate(domain_, boundary_,
                           atomStorage_, bondStorage_,
                           #ifdef INTER_ANGLE
                           angleStorage_,
                           #endif
                           #ifdef INTER_DIHEDRAL
                           dihedralStorage_,
                           #endif
                           buffer_);
      #endif
      exchanger_.associate(domain_, boundary_, atomStorage_, buffer_);
      exchanger_.addGroupExchanger(bondStorage_);
      #ifdef INTER_ANGLE
      exchanger_.addGroupExchanger(angleStorage_);
      #endif
      #ifdef INTER_DIHEDRAL
      exchanger_.addGroupExchanger(dihedralStorage_);
      #endif

      fileMasterPtr_ = new FileMaster;
      energyEnsemblePtr_  = new EnergyEnsemble;
      boundaryEnsemblePtr_ = new BoundaryEnsemble;

      diagnosticManagerPtr_ = new DiagnosticManager(*this);
   }

   /*
   * Destructor.
   */
   Simulation::~Simulation()
   {
      if (pairPotentialPtr_) {
         delete pairPotentialPtr_;
      }
      if (bondPotentialPtr_) {
         delete bondPotentialPtr_;
      }
      #ifdef INTER_ANGLE
      if (anglePotentialPtr_) {
         delete anglePotentialPtr_;
      }
      #endif
      #ifdef INTER_DIHEDRAL
      if (dihedralPotentialPtr_) {
         delete dihedralPotentialPtr_;
      }
      #endif
      #ifdef INTER_EXTERNAL
      if (externalPotentialPtr_) {
         delete externalPotentialPtr_;
      }
      #endif
      if (energyEnsemblePtr_) {
         delete energyEnsemblePtr_;
      }
      if (boundaryEnsemblePtr_) {
         delete boundaryEnsemblePtr_;
      }
      if (integratorPtr_) {
         delete integratorPtr_;
      }
      if (fileMasterPtr_) {
         delete fileMasterPtr_;
      }
      if (configIoPtr_) {
         delete configIoPtr_;
      }
      if (serializeConfigIoPtr_) {
         delete serializeConfigIoPtr_;
      }
      if (diagnosticManagerPtr_) {
         delete diagnosticManagerPtr_;
      }
      #ifndef DDMD_NOPAIR
      if (pairFactoryPtr_) {
         delete pairFactoryPtr_;
      }
      #endif
      if (bondFactoryPtr_) {
         delete bondFactoryPtr_;
      }
      #ifdef INTER_ANGLE
      if (angleFactoryPtr_) {
         delete angleFactoryPtr_;
      }
      #endif
      #ifdef INTER_DIHEDRAL
      if (dihedralFactoryPtr_) {
         delete dihedralFactoryPtr_;
      }
      #endif
      #ifdef INTER_EXTERNAL
      if (externalFactoryPtr_) {
         delete externalFactoryPtr_;
      }
      #endif
      if (integratorFactoryPtr_) {
         delete integratorFactoryPtr_;
      }
      if (configIoFactoryPtr_) {
         delete configIoFactoryPtr_;
      }

      #ifdef UTIL_MPI
      if (logFile_.is_open()) logFile_.close();
      #endif
   }

   /*
   * Process command line options.
   */
   void Simulation::setOptions(int argc, char * const * argv)
   {
      bool  eFlag = false;
      bool  rFlag = false;
      bool  sFlag = false;
      char* sArg;
      char* rArg;
      int  nSystem = 1;

      // Read command-line arguments
      int c;
      opterr = 0;
      while ((c = getopt(argc, argv, "er:s:")) != -1) {
         switch (c) {
         case 'e':
           eFlag = true;
           break;
         case 'r':
           rFlag = true;
           rArg  = optarg;
           break;
         case 's':
           sFlag = true;
           sArg  = optarg;
           nSystem = atoi(sArg);
           break;
         case '?':
           std::cout << "Unknown option -" << optopt << std::endl;
         }
      }

      // If option -s, split the communicator
      if (sFlag) {
         if (nSystem <= 1) {
            UTIL_THROW("nSystem must be greater than 1");
         }
         int worldRank = communicator_.Get_rank();
         int worldSize = communicator_.Get_size();
         if (worldSize % nSystem != 0) {
            UTIL_THROW("World communicator size not a multiple of nSystem");
         }

         // Split the communicator
         int systemSize = worldSize/nSystem;
         int systemId  = worldRank/systemSize;
         communicator_ = communicator_.Split(systemId, worldRank);

         // Set param and grid communicators
         setIoCommunicator(communicator_);
         domain_.setGridCommunicator(communicator_);
         fileMaster().setDirectoryId(systemId);

         // Set log file for processor n to a new file named "n/log".
         // Relies on initialization of outputPrefix to empty string "".
         fileMaster().openOutputFile("log", logFile_);
         Log::setFile(logFile_);
      }

      // If option -e, enable echoing of parameters as they are read
      if (eFlag) {
         ParamComponent::setEcho(true);
      }

      // If option -r, load state from a restart file.
      if (rFlag) {
         if (isIoProcessor()) {
            Log::file() << "Begin reading restart, base file name "
                        << std::string(rArg) << std::endl;
         }
         load(std::string(rArg));
         if (isIoProcessor()) {
            Log::file() << std::endl;
         }
      }

   }

   /*
   *  Read parameters from default parameter file.
   */
   void Simulation::readParam()
   {  readParam(fileMaster().paramFile()); }

   /*
   * Read parameter block, including begin and end.
   */
   void Simulation::readParam(std::istream& in)
   {
      if (isRestarting_) {
         if (isInitialized_) {
            return;
         }
      }
      readBegin(in, className().c_str());
      readParameters(in);
      readEnd(in);
   }

   /*
   * Read parameters, allocate memory and initialize.
   */
   void Simulation::readParameters(std::istream& in)
   {
      // Read Domain and FileMaster
      readParamComposite(in, domain_);
      readFileMaster(in);

      // Read numbers of types
      read<int>(in, "nAtomType", nAtomType_);
      read<int>(in, "nBondType", nBondType_);
      #ifdef INTER_ANGLE
      read<int>(in, "nAngleType", nAngleType_);
      #endif
      #ifdef INTER_DIHEDRAL
      read<int>(in, "nDihedralType", nDihedralType_);
      #endif
      #ifdef INTER_EXTERNAL
      read<bool>(in, "hasExternal", hasExternal_);
      #endif

      // Read array of atom type descriptors
      atomTypes_.allocate(nAtomType_);
      for (int i = 0; i < nAtomType_; ++i) {
         atomTypes_[i].setId(i);
      }
      readDArray<AtomType>(in, "atomTypes", atomTypes_, nAtomType_);

      // Read storage capacities
      readParamComposite(in, atomStorage_);
      readParamComposite(in, bondStorage_);
      #ifdef INTER_ANGLE
      if (nAngleType_) {
         readParamComposite(in, angleStorage_);
      }
      #endif
      #ifdef INTER_DIHEDRAL
      if (nDihedralType_) {
         readParamComposite(in, dihedralStorage_);
      }
      #endif

      readParamComposite(in, buffer_);
      readPotentialStyles(in);

      // Create and read potential energy classes

      #ifndef DDMD_NOPAIR
      // Pair Potential
      assert(pairPotentialPtr_ == 0);
      pairPotentialPtr_ = pairFactory().factory(pairStyle());
      if (!pairPotentialPtr_) {
         UTIL_THROW("Unknown pairStyle");
      }
      pairPotentialPtr_->setNAtomType(nAtomType_);
      readParamComposite(in, *pairPotentialPtr_);
      pairPotentialPtr_->setReverseUpdateFlag(reverseUpdateFlag_);
      #endif

      // Bond Potential
      assert(bondPotentialPtr_ == 0);
      bondPotentialPtr_ = bondFactory().factory(bondStyle());
      if (!bondPotentialPtr_) {
         UTIL_THROW("Unknown bondStyle");
      }
      bondPotentialPtr_->setNBondType(nBondType_);
      readParamComposite(in, *bondPotentialPtr_);

      #ifdef INTER_ANGLE
      // Angle potential
      if (nAngleType_) {
         assert(anglePotentialPtr_ == 0);
         anglePotentialPtr_ = angleFactory().factory(angleStyle());
         if (!anglePotentialPtr_) {
            UTIL_THROW("Unknown angleStyle");
         }
         anglePotentialPtr_->setNAngleType(nAngleType_);
         readParamComposite(in, *anglePotentialPtr_);
      }
      #endif

      #ifdef INTER_DIHEDRAL
      // Dihedral potential
      if (nDihedralType_) {
         assert(dihedralPotentialPtr_ == 0);
         dihedralPotentialPtr_ = dihedralFactory().factory(dihedralStyle());
         if (!dihedralPotentialPtr_) {
            UTIL_THROW("Unknown dihedralStyle");
         }
         dihedralPotentialPtr_->setNDihedralType(nDihedralType_);
         readParamComposite(in, *dihedralPotentialPtr_);
      }
      #endif

      #ifdef INTER_EXTERNAL
      // External potential
      if (hasExternal_) {
         assert(externalPotentialPtr_ == 0);
         externalPotentialPtr_ = externalFactory().factory(externalStyle());
         externalPotentialPtr_->setNAtomType(nAtomType_);
         readParamComposite(in, *externalPotentialPtr_);
      }
      #endif

      readEnsembles(in);

      // Integrator
      std::string className;
      bool isEnd;
      assert(integratorPtr_ == 0);
      integratorPtr_ =
         integratorFactory().readObject(in, *this, className, isEnd);
      if (!integratorPtr_) {
         std::string msg("Unknown Integrator subclass name: ");
         msg += className;
         UTIL_THROW("msg.c_str()");
      }

      readParamComposite(in, random_);
      readParamComposite(in, *diagnosticManagerPtr_);

      // Finished reading parameter file. Now finish initialization:

      exchanger_.setPairCutoff(pairPotentialPtr_->cutoff());
      exchanger_.allocate();

      // Set signal observers (i.e., call-back functions for Signal::notify)
      modifySignal().addObserver(*this, &Simulation::unsetKineticEnergy);
      modifySignal().addObserver(*this, &Simulation::unsetKineticStress);
      modifySignal().addObserver(*this, &Simulation::unsetPotentialEnergies);
      modifySignal().addObserver(*this, &Simulation::unsetVirialStress);

      velocitySignal().addObserver(*this, &Simulation::unsetKineticEnergy);
      velocitySignal().addObserver(*this, &Simulation::unsetKineticStress);

      positionSignal().addObserver(*this, &Simulation::unsetPotentialEnergies);
      positionSignal().addObserver(*this, &Simulation::unsetVirialStress);
      if (nBondType_) {
         void (BondStorage::*memberPtr)() = &BondStorage::unsetNTotal;
         exchangeSignal().addObserver(bondStorage_, memberPtr);
      }
      #ifdef INTER_ANGLE
      if (nAngleType_) {
         void (AngleStorage::*memberPtr)() = &AngleStorage::unsetNTotal;
         exchangeSignal().addObserver(angleStorage_, memberPtr);
      }
      #endif
      #ifdef INTER_DIHEDRAL
      if (nDihedralType_) {
         void (DihedralStorage::*memberPtr)() = &DihedralStorage::unsetNTotal;
         exchangeSignal().addObserver(dihedralStorage_, memberPtr);
      }
      #endif

      isInitialized_ = true;
   }

   /*
   * Read parameters, allocate memory and initialize.
   */
   void Simulation::loadParameters(Serializable::IArchive& ar)
   {
      if (isInitialized_) {
         UTIL_THROW("Error: Called loadParameters when already initialized");
      }
      if (!hasIoCommunicator()) {
         UTIL_THROW("Error: No IoCommunicator is set");
      }
      isRestarting_ = true;

      loadParamComposite(ar, domain_);
      loadFileMaster(ar);

      // Load types
      loadParameter<int>(ar, "nAtomType", nAtomType_);
      loadParameter<int>(ar, "nBondType", nBondType_);
      #ifdef INTER_ANGLE
      loadParameter<int>(ar, "nAngleType", nAngleType_);
      #endif
      #ifdef INTER_DIHEDRAL
      loadParameter<int>(ar, "nDihedralType", nDihedralType_);
      #endif
      #ifdef INTER_EXTERNAL
      loadParameter<bool>(ar, "hasExternal", hasExternal_);
      #endif
      atomTypes_.allocate(nAtomType_);
      for (int i = 0; i < nAtomType_; ++i) {
         atomTypes_[i].setId(i);
      }
      loadDArray<AtomType>(ar, "atomTypes", atomTypes_, nAtomType_);

      // Load storage capacities
      loadParamComposite(ar, atomStorage_);
      loadParamComposite(ar, bondStorage_);
      #ifdef INTER_ANGLE
      if (nAngleType_) {
         loadParamComposite(ar, angleStorage_);
      }
      #endif
      #ifdef INTER_DIHEDRAL
      if (nDihedralType_) {
         loadParamComposite(ar, dihedralStorage_);
      }
      #endif
      loadParamComposite(ar, buffer_);

      // Load potentials styles and parameters
      loadPotentialStyles(ar);
      #ifndef DDMD_NOPAIR
      // Pair Potential
      assert(pairPotentialPtr_ == 0);
      pairPotentialPtr_ = pairFactory().factory(pairStyle());
      if (!pairPotentialPtr_) {
         UTIL_THROW("Unknown pairStyle");
      }
      pairPotentialPtr_->setNAtomType(nAtomType_);
      loadParamComposite(ar, *pairPotentialPtr_);
      pairPotentialPtr_->setReverseUpdateFlag(reverseUpdateFlag_);
      #endif
      // Bond Potential
      assert(bondPotentialPtr_ == 0);
      bondPotentialPtr_ = bondFactory().factory(bondStyle());
      if (!bondPotentialPtr_) {
         UTIL_THROW("Unknown bondStyle");
      }
      bondPotentialPtr_->setNBondType(nBondType_);
      loadParamComposite(ar, *bondPotentialPtr_);
      #ifdef INTER_ANGLE
      // Angle potential
      assert(anglePotentialPtr_ == 0);
      if (nAngleType_) {
         anglePotentialPtr_ = angleFactory().factory(angleStyle());
         if (!anglePotentialPtr_) {
            UTIL_THROW("Unknown angleStyle");
         }
         anglePotentialPtr_->setNAngleType(nAngleType_);
         loadParamComposite(ar, *anglePotentialPtr_);
      }
      #endif
      #ifdef INTER_DIHEDRAL
      // Dihedral potential
      assert(dihedralPotentialPtr_ == 0);
      if (nDihedralType_) {
         dihedralPotentialPtr_ = dihedralFactory().factory(dihedralStyle());
         if (!dihedralPotentialPtr_) {
            UTIL_THROW("Unknown dihedralStyle");
         }
         dihedralPotentialPtr_->setNDihedralType(nDihedralType_);
         loadParamComposite(ar, *dihedralPotentialPtr_);
      }
      #endif
      #ifdef INTER_EXTERNAL
      // External potential
      if (hasExternal_) {
         externalPotentialPtr_ = externalFactory().factory(externalStyle());
         externalPotentialPtr_->setNAtomType(nAtomType_);
         loadParamComposite(ar, *externalPotentialPtr_);
      }
      #endif

      loadEnsembles(ar);

      // Integrator
      assert(integratorPtr_ == 0);
      std::string className;
      bool isEnd;
      integratorPtr_ =
         integratorFactory().loadObject(ar, *this, className);
      if (!integratorPtr_) {
         std::string msg("Unknown Integrator subclass name: ");
         msg += className;
         UTIL_THROW("msg.c_str()");
      }

      loadParamComposite(ar, random_);
      loadParamComposite(ar, *diagnosticManagerPtr_);

      // Finished loading data from archive. Finish initialization:

      exchanger_.setPairCutoff(pairPotentialPtr_->cutoff());
      exchanger_.allocate();

      // Set signal observers (i.e., call-back functions for Signal::notify)
      modifySignal().addObserver(*this, &Simulation::unsetKineticEnergy);
      modifySignal().addObserver(*this, &Simulation::unsetKineticStress);
      modifySignal().addObserver(*this, &Simulation::unsetPotentialEnergies);
      modifySignal().addObserver(*this, &Simulation::unsetVirialStress);

      velocitySignal().addObserver(*this, &Simulation::unsetKineticEnergy);
      velocitySignal().addObserver(*this, &Simulation::unsetKineticStress);

      positionSignal().addObserver(*this, &Simulation::unsetPotentialEnergies);
      positionSignal().addObserver(*this, &Simulation::unsetVirialStress);
      if (nBondType_) {
         void (BondStorage::*memberPtr)() = &BondStorage::unsetNTotal;
         exchangeSignal().addObserver(bondStorage_, memberPtr);
      }
      #ifdef INTER_ANGLE
      if (nAngleType_) {
         void (AngleStorage::*memberPtr)() = &AngleStorage::unsetNTotal;
         exchangeSignal().addObserver(angleStorage_, memberPtr);
      }
      #endif
      #ifdef INTER_DIHEDRAL
      if (nDihedralType_) {
         void (DihedralStorage::*memberPtr)() = &DihedralStorage::unsetNTotal;
         exchangeSignal().addObserver(dihedralStorage_, memberPtr);
      }
      #endif

      isInitialized_ = true;

      // Load the configuration (boundary + positions + groups)
      serializeConfigIo().loadConfig(ar, maskedPairPolicy_);

      // There are no ghosts yet, so exchange.
      exchanger_.exchange();
      isValid();
   }

   // ---- Serialization -----------------------------------------------

   /*
   * Load state from a restart file (open file, call load, close file)
   */
   void Simulation::load(const std::string& filename)
   {
      if (isInitialized_) {
         UTIL_THROW("Error: Called load when already initialized");
      }
      if (!hasIoCommunicator()) {
         UTIL_THROW("Error: No IoCommunicator is set");
      }

      Serializable::IArchive ar;
      if (isIoProcessor()) {
         fileMaster().openRestartIFile(filename, ".rst", ar.file());
      }
      // Call ParamComposite::load(), which calls Simulation::loadParameters()
      load(ar); 
      if (isIoProcessor()) {
         ar.file().close();
      }

      // Set default command (*.cmd) file name = filename + ".cmd".
      // This will be used by Simulation::readCommands().
      std::string commandFileName = filename + ".cmd";
      fileMaster().setCommandFileName(commandFileName);
   }

   /*
   * Serialize internal state to an archive.
   *
   * Call only on IoProcessor of IoCommunicator.
   */
   void Simulation::save(Serializable::OArchive& ar)
   {
      domain_.save(ar);
      saveFileMaster(ar);

      // Read types
      ar << nAtomType_;
      ar << nBondType_;
      #ifdef INTER_ANGLE
      ar << nAngleType_;
      #endif
      #ifdef INTER_DIHEDRAL
      ar << nDihedralType_;
      #endif
      #ifdef INTER_EXTERNAL
      ar << hasExternal_;
      #endif
      ar << atomTypes_;

      // Read storage capacities
      atomStorage_.save(ar);
      bondStorage_.save(ar);
      #ifdef INTER_ANGLE
      if (nAngleType_) {
         angleStorage_.save(ar);
      }
      #endif
      #ifdef INTER_DIHEDRAL
      if (nDihedralType_) {
         dihedralStorage_.save(ar);
      }
      #endif

      buffer_.save(ar);

      // Potentials energy styles and parameters
      savePotentialStyles(ar);
      #ifndef DDMD_NOPAIR
      assert(pairPotentialPtr_);
      pairPotentialPtr_->save(ar);
      #endif
      if (bondPotentialPtr_) {
         assert(bondPotentialPtr_);
         bondPotentialPtr_->save(ar);
      }
      #ifdef INTER_ANGLE
      if (nAngleType_) {
         assert(anglePotentialPtr_);
         anglePotentialPtr_->save(ar);
      }
      #endif
      #ifdef INTER_DIHEDRAL
      if (nDihedralType_) {
         assert(dihedralPotentialPtr_);
         dihedralPotentialPtr_->save(ar);
      }
      #endif
      #ifdef INTER_EXTERNAL
      // External potential
      if (hasExternal_) {
         assert(externalPotentialPtr_);
         externalPotentialPtr_->save(ar);
      }
      #endif

      // Save ensembles and integrator
      saveEnsembles(ar);
      std::string name = integrator().className();
      ar << name;
      integrator().save(ar);

      random_.save(ar);
      diagnosticManagerPtr_->save(ar);
   }

   /*
   * Save state to file (open file, call save(), close file).
   */
   void Simulation::save(const std::string& filename)
   {
      // Update statistics (call on all processors).
      atomStorage_.computeStatistics(domain_.communicator());
      bondStorage_.computeStatistics(domain_.communicator());
      #ifdef INTER_ANGLE
      angleStorage_.computeStatistics(domain_.communicator());
      #endif
      #ifdef INTER_DIHEDRAL
      dihedralStorage_.computeStatistics(domain_.communicator());
      #endif
      buffer_.computeStatistics(domain_.communicator());
      if (integratorPtr_) {
         integrator().computeStatistics();
      }

      // Save parameters (only on ioProcessor)
      Serializable::OArchive ar;
      if (isIoProcessor()) {
         fileMaster().openRestartOFile(filename, ".rst", ar.file());
         save(ar);
      }

      // Save configuration (call on all processors)
      serializeConfigIo().saveConfig(ar);

      if (isIoProcessor()) {
         ar.file().close();
      }
   }

   // --- Protected read, load and save methods ------------------------

   /*
   * Read the FileMaster parameters.
   */
   void Simulation::readFileMaster(std::istream &in)
   {  readParamComposite(in, *fileMasterPtr_); }

   /*
   * If no FileMaster exists, create and load one.
   */
   void Simulation::loadFileMaster(Serializable::IArchive& ar)
   {  loadParamComposite(ar, *fileMasterPtr_); }

   /*
   * If createdFileMaster_, save to archive.
   *
   * A Simulation normally creates its own FileMaster only in unit tests.
   * In a simulation, the FileMaster is normally set to that of either
   * a parent Simulation or that of another Simulation from which this was
   * copied.
   */
   void Simulation::saveFileMaster(Serializable::OArchive& ar)
   {  fileMasterPtr_->save(ar); }

   /*
   * Read potential style strings and maskedPairPolicy.
   */
   void Simulation::readPotentialStyles(std::istream &in)
   {
      #ifndef DDMD_NOPAIR
      read<std::string>(in, "pairStyle", pairStyle_);
      #endif
      if (nBondType_) {
         read<std::string>(in, "bondStyle", bondStyle_);
      }
      #ifdef INTER_ANGLE
      if (nAngleType_) {
         read<std::string>(in, "angleStyle", angleStyle_);
      }
      #endif
      #ifdef INTER_DIHEDRAL
      if (nDihedralType_) {
         read<std::string>(in, "dihedralStyle", dihedralStyle_);
      }
      #endif
      #ifdef INTER_EXTERNAL
      if (hasExternal_) {
         read<std::string>(in, "externalStyle", externalStyle_);
      }
      #endif

      // Read policy regarding whether to excluded pair interactions
      // between covalently bonded pairs.
      read<MaskPolicy>(in, "maskedPairPolicy", maskedPairPolicy_);

      // Reverse communication (true) or not (false)?
      read<bool>(in, "reverseUpdateFlag", reverseUpdateFlag_);

   }

   /*
   * Load potential style strings.
   */
   void Simulation::loadPotentialStyles(Serializable::IArchive& ar)
   {
      #ifndef INTER_NOPAIR
      loadParameter<std::string>(ar, "pairStyle", pairStyle_);
      #endif
      if (nBondType()) {
         loadParameter<std::string>(ar, "bondStyle", bondStyle_);
      }
      #ifdef INTER_ANGLE
      if (nAngleType()) {
         loadParameter<std::string>(ar, "angleStyle", angleStyle_);
      }
      #endif
      #ifdef INTER_DIHEDRAL
      if (nDihedralType()) {
         loadParameter<std::string>(ar, "dihedralStyle", dihedralStyle_);
      }
      #endif
      #ifdef INTER_EXTERNAL
      if (hasExternal()) {
         loadParameter<std::string>(ar, "externalStyle", externalStyle_);
      }
      #endif
      loadParameter<MaskPolicy>(ar, "maskedPairPolicy", maskedPairPolicy_);
      loadParameter<bool>(ar, "reverseUpdateFlag", reverseUpdateFlag_);

      isInitialized_ = true;
   }

   /*
   * Save potential style strings.
   */
   void Simulation::savePotentialStyles(Serializable::OArchive& ar)
   {
      #ifndef INTER_NOPAIR
      ar << pairStyle_;
      #endif
      if (nBondType()) {
         ar << bondStyle_;
      }
      #ifdef INTER_ANGLE
      if (nAngleType()) {
         ar << angleStyle_;
      }
      #endif
      #ifdef INTER_DIHEDRAL
      if (nDihedralType()) {
         ar << dihedralStyle_;
      }
      #endif
      #ifdef INTER_EXTERNAL
      if (hasExternal()) {
         ar << externalStyle_;
      }
      #endif
      ar << maskedPairPolicy_;
      ar << reverseUpdateFlag_;
   }

   /*
   * Read EnergyEnsemble and BoundaryEnsemble
   */
   void Simulation::readEnsembles(std::istream &in)
   {
      readParamComposite(in, *energyEnsemblePtr_);
      readParamComposite(in, *boundaryEnsemblePtr_);
   }

   /*
   * Load EnergyEnsemble and BoundaryEnsemble.
   */
   void Simulation::loadEnsembles(Serializable::IArchive& ar)
   {
      loadParamComposite(ar, *energyEnsemblePtr_);
      loadParamComposite(ar, *boundaryEnsemblePtr_);
   }

   /*
   * Save EnergyEnsemble and BoundaryEnsemble.
   */
   void Simulation::saveEnsembles(Serializable::OArchive& ar)
   {
      energyEnsemblePtr_->save(ar);
      boundaryEnsemblePtr_->save(ar);
   }

   // --- readCommands and run-time actions  ---------------------------

   /*
   * Read and implement commands in an input script.
   */
   void Simulation::readCommands(std::istream &in)
   {
      std::string  command;
      std::string  filename;
      std::ifstream  inputFile;
      std::ofstream  outputFile;

      #ifdef UTIL_MPI
      std::stringstream  inBuffer;
      std::string  line;
      #else
      std::istream&  inBuffer = in;
      #endif

      bool readNext = true;
      while (readNext) {

         // Precondition: Check atomic coordinate system
         if (atomStorage().isCartesian()) {
            UTIL_THROW("Error: Storage set for Cartesian atom coordinates");
         }

         #ifdef UTIL_MPI
         // Read and broadcast command line
         if (!hasIoCommunicator() || isIoProcessor()) {
            getNextLine(in, line);
            Log::file() << line << std::endl;
         }
         if (hasIoCommunicator()) {
            bcast<std::string>(domain_.communicator(), line, 0);
         }
         #else // not UTIL_MPI
         getNextLine(in, line);
         Log::file() << line << std::endl;
         #endif // not UTIL_MPI
         inBuffer.clear();
         for (unsigned i=0; i < line.size(); ++i) {
            inBuffer.put(line[i]);
         }

         inBuffer >> command;
         if (isRestarting_) {

            if (command == "RESTART") {
               int endStep;
               inBuffer >> endStep;
               int iStep = integrator().iStep();
               int nStep = endStep - iStep;
               if (isIoProcessor()) {
                  Log::file() << "Running from  iStep =" << iStep << " to "
                              << endStep << std::endl;
               }
               integrator().run(nStep);
               isRestarting_ = false;
            } else {
               UTIL_THROW("Missing RESTART command when restarting");
            }

         } else {

            if (command == "READ_CONFIG") {
               // Read configuration (boundary + positions + topology).
               // Use current ConfigIo, if set, or create default.
               inBuffer >> filename;
               readConfig(filename);
            } else
            if (command == "THERMALIZE") {
               double temperature;
               inBuffer >> temperature;
               setBoltzmannVelocities(temperature);
            } else
            if (command == "SIMULATE") {
               int nStep;
               inBuffer >> nStep;
               if (domain_.isMaster()) {
                  Log::file() << std::endl;
               }
               integrator().run(nStep);
            } else
            if (command == "OUTPUT_DIAGNOSTICS") {
               diagnosticManager().output();
            } else
            if (command == "OUTPUT_INTEGRATOR_STATS") {
               // Output statistics about time usage during simulation.
               integrator().computeStatistics();
               if (domain_.isMaster()) {
                  integrator().outputStatistics(Log::file());
               }
            } else
            if (command == "OUTPUT_EXCHANGER_STATS") {
               // Output detailed statistics about time usage by Exchanger.
               integrator().computeStatistics();
               #ifdef UTIL_MPI
               exchanger().timer().reduce(domain().communicator());
               #endif
               if (domain_.isMaster()) {
                  int iStep = integrator().iStep();
                  double time  = integrator().time();
                  exchanger_.outputStatistics(Log::file(), time, iStep);
               }
            } else
            if (command == "OUTPUT_MEMORY_STATS") {
               // Output statistics about memory usage during simulation.
               // Also clears statistics after printing output
               atomStorage().computeStatistics(domain_.communicator());
               bondStorage().computeStatistics(domain_.communicator());
               #ifdef INTER_ANGLE
               angleStorage().computeStatistics(domain_.communicator());
               #endif
               #ifdef INTER_DIHEDRAL
               dihedralStorage().computeStatistics(domain_.communicator());
               #endif
               buffer().computeStatistics(domain_.communicator());
               pairPotential().pairList()
                              .computeStatistics(domain_.communicator());
               if (domain_.isMaster()) {
                  atomStorage().outputStatistics(Log::file());
                  bondStorage().outputStatistics(Log::file());
                  buffer().outputStatistics(Log::file());
                  pairPotential().pairList().outputStatistics(Log::file());
                  Log::file() << std::endl;
               }

               atomStorage().clearStatistics();
               bondStorage().clearStatistics();
               #ifdef INTER_ANGLE
               angleStorage().clearStatistics();
               #endif
               #ifdef INTER_DIHEDRAL
               dihedralStorage().clearStatistics();
               #endif
               buffer().clearStatistics();
               pairPotential().pairList().clearStatistics();

            } else
            if (command == "CLEAR_INTEGRATOR") {
               // Clear timing, memory statistics, diagnostic accumulators.
               // Also resets integrator iStep() to zero
               integrator().clear();
            } else
            if (command == "WRITE_CONFIG") {
               // Write current configuration to file.
               // Use current ConfigIo, if set.
               // Otherwise, create instance of default ConfigIo class.
               inBuffer >> filename;
               writeConfig(filename);
            } else
            if (command == "WRITE_PARAM") {
               // Write params file, using current variable values.
               inBuffer >> filename;
               fileMaster().openOutputFile(filename, outputFile);
               writeParam(outputFile);
               outputFile.close();
            } else
            if (command == "SET_CONFIG_IO") {
               // Create a ConfigIo object of specified subclass.
               // This gives file format for later reads and writes.
               std::string classname;
               inBuffer >> classname;
               setConfigIo(classname);
            } else
            if (command == "SET_PAIR") {
               // Modify one parameter of a pair interaction.
               std::string paramName;
               int typeId1, typeId2;
               double value;
               inBuffer >> paramName >> typeId1 >> typeId2 >> value;
               pairPotential().set(paramName, typeId1, typeId2, value);
            } else
            if (command == "SET_BOND") {
               // Modify one parameter of a bond interaction.
               std::string paramName;
               int typeId;
               double value;
               inBuffer >> paramName >> typeId >> value;
               bondPotential().set(paramName, typeId, value);
            } else
            #ifdef INTER_ANGLE
            if (command == "SET_ANGLE" && nAngleType_) {
               // Modify one parameter of an angle interaction.
               std::string paramName;
               int typeId;
               double value;
               inBuffer >> paramName >> typeId >> value;
               anglePotential().set(paramName, typeId, value);
            } else
            #endif
            #ifdef INTER_DIHEDRAL
            if (command == "SET_DIHEDRAL" && nDihedralType_) {
               // Modify one parameter of a dihedral interaction.
               std::string paramName;
               int typeId;
               double value;
               inBuffer >> paramName >> typeId >> value;
               dihedralPotential().set(paramName, typeId, value);
            } else
            #endif
            if (command == "FINISH") {
               // Terminate loop over commands.
               readNext = false;
            } else {
               Log::file() << "Error: Unknown command  " << std::endl;
               readNext = false;
            }
         }

      }
   }

   /*
   * Read and implement commands from the default command file.
   */
   void Simulation::readCommands()
   {  readCommands(fileMaster().commandFile()); }

   /*
   * Set flag to specify if reverse communication is enabled.
   */
   void Simulation::setReverseUpdateFlag(bool reverseUpdateFlag)
   {
      reverseUpdateFlag_ = reverseUpdateFlag;
      if (pairPotentialPtr_) {
         pairPotentialPtr_->setReverseUpdateFlag(reverseUpdateFlag);
      }
   }

   /*
   * Choose velocities from a Boltzmann distribution.
   */
   void Simulation::setBoltzmannVelocities(double temperature)
   {
      double scale = sqrt(temperature);
      AtomIterator atomIter;
      int i;

      atomStorage_.begin(atomIter);
      for( ; atomIter.notEnd(); ++atomIter){
         for (i = 0; i < Dimension; ++i) {
            atomIter->velocity()[i] = scale*random_.gaussian();
         }
      }
      velocitySignal().notify();
   }

   /*
   * Set forces on all local atoms to zero.
   * If reverseUpdateFlag(), also zero ghost atom forces.
   */
   void Simulation::zeroForces()
   {
      // Zero local atoms
      AtomIterator atomIter;
      atomStorage_.begin(atomIter);
      for( ; atomIter.notEnd(); ++atomIter){
         atomIter->force().zero();
      }

      // If using reverse communication, zero ghost atoms
      if (reverseUpdateFlag_) {
         GhostIterator ghostIter;
         atomStorage_.begin(ghostIter);
         for( ; ghostIter.notEnd(); ++ghostIter){
            ghostIter->force().zero();
         }
      }
   }

   /*
   * Compute forces for all atoms.
   */
   void Simulation::computeForces()
   {
      zeroForces();
      pairPotential().computeForces();
      bondPotential().computeForces();
      #ifdef INTER_ANGLE
      if (nAngleType_) {
         anglePotential().computeForces();
      }
      #endif
      #ifdef INTER_DIHEDRAL
      if (nDihedralType_) {
         dihedralPotential().computeForces();
      }
      #endif
      #ifdef INTER_EXTERNAL
      if (hasExternal_) {
         externalPotential().computeForces();
      }
      #endif

      // Reverse communication (if any)
      if (reverseUpdateFlag_) {
         exchanger_.reverseUpdate();
      }
   }

   /*
   * Compute forces for all atoms and virial stress contributions.
   */
   void Simulation::computeForcesAndVirial()
   {
      zeroForces();
      pairPotential().computeForcesAndStress(domain_.communicator());
      bondPotential().computeForcesAndStress(domain_.communicator());
      #ifdef INTER_ANGLE
      if (nAngleType_) {
         anglePotential().computeForcesAndStress(domain_.communicator());
      }
      #endif
      #ifdef INTER_DIHEDRAL
      if (nDihedralType_) {
         dihedralPotential().computeForcesAndStress(domain_.communicator());
      }
      #endif
      #ifdef INTER_EXTERNAL
      if (hasExternal_) {
         externalPotential().computeForces();
      }
      #endif

      // Reverse communication (if any)
      if (reverseUpdateFlag_) {
         exchanger_.reverseUpdate();
      }
   }

   // --- Kinetic Energy methods ---------------------------------------

   /*
   * Calculate total kinetic energy (call on all processors).
   */
   void Simulation::computeKineticEnergy()
   {
      // Do nothing if kinetic energy is already set.
      if (kineticEnergy_.isSet()) return;

      double localEnergy = 0.0;
      double mass;
      int typeId;

      // Add kinetic energies of local atoms on this processor
      AtomIterator atomIter;
      atomStorage_.begin(atomIter);
      for( ; atomIter.notEnd(); ++atomIter){
         typeId = atomIter->typeId();
         mass   = atomTypes_[typeId].mass();
         localEnergy += mass*(atomIter->velocity().square());
      }
      localEnergy = 0.5*localEnergy;

      #ifdef UTIL_MPI
      // Sum values from all processors.
      double totalEnergy = 0.0;
      domain_.communicator().Reduce(&localEnergy, &totalEnergy, 1,
                                    MPI::DOUBLE, MPI::SUM, 0);
      if (domain_.communicator().Get_rank() != 0) {
         totalEnergy = 0.0;
      }
      kineticEnergy_.set(totalEnergy);
      #else
      kineticEnergy_.set(localEnergy);
      #endif
   }

   /*
   * Return total kinetic energy (on master processor).
   *
   * Call only on master processor.
   */
   double Simulation::kineticEnergy()
   {  return kineticEnergy_.value(); }

   /*
   * Mark kinetic energy as unknown.
   */
   void Simulation::unsetKineticEnergy()
   {  kineticEnergy_.unset(); }

   // --- Kinetic Stress methods ---------------------------------------

   /*
   * Compute total kinetic stress, store on master proc.
   *
   * Call on all processors.
   */
   void Simulation::computeKineticStress()
   {

      // Do nothing and return if kinetic stress is already set
      if (kineticStress_.isSet()) return;

      Tensor localStress;
      double  mass;
      Vector  velocity;
      int typeId, i, j;

      localStress.zero();

      // For each local atoms on this processor, stress(i, j) += m*v[i]*v[j]
      AtomIterator atomIter;
      atomStorage_.begin(atomIter);
      for( ; atomIter.notEnd(); ++atomIter){
         typeId = atomIter->typeId();
         velocity = atomIter->velocity();
         mass = atomTypes_[typeId].mass();
         for (i = 0; i < Dimension; ++i) {
            for (j = 0; j < Dimension; ++j) {
               localStress(i, j) += mass*velocity[i]*velocity[j];
            }
         }
      }

      // Divide dyad sum by system volume
      localStress /= boundary().volume();

      #ifdef UTIL_MPI
      // Sum values from all processors
      Tensor totalStress;
      domain_.communicator().Reduce(&localStress(0, 0), &totalStress(0, 0),
                                    Dimension*Dimension, MPI::DOUBLE, MPI::SUM, 0);
      if (domain_.communicator().Get_rank() != 0) {
         totalStress.zero();
      }
      kineticStress_.set(totalStress);
      #else
      kineticStress_.set(localStress);
      #endif
   }

   /*
   * Return total kinetic stress (on master processor).
   */
   Tensor Simulation::kineticStress() const
   {  return kineticStress_.value(); }

   /*
   * Return total kinetic pressure (on master processor).
   */
   double Simulation::kineticPressure() const
   {
      double pressure = 0.0;
      for (int i = 0; i < Dimension; ++i) {
         pressure += kineticStress_.value()(i, i);
      }
      pressure /= 3.0;
      return pressure;
   }

   /*
   * Mark kinetic stress as unknown.
   */
   void Simulation::unsetKineticStress()
   {  kineticStress_.unset(); }

   // --- Potential Energy Methods -------------------------------------

   #ifdef UTIL_MPI
   /*
   * Compute all potential energy contributions.
   */
   void Simulation::computePotentialEnergies()
   {
      pairPotential().computeEnergy(domain_.communicator());
      bondPotential().computeEnergy(domain_.communicator());
      #ifdef INTER_ANGLE
      if (nAngleType_) {
         anglePotential().computeEnergy(domain_.communicator());
      }
      #endif
      #ifdef INTER_DIHEDRAL
      if (nDihedralType_) {
         dihedralPotential().computeEnergy(domain_.communicator());
      }
      #endif
      #ifdef INTER_EXTERNAL
      if (hasExternal_) {
         externalPotential().computeEnergy(domain_.communicator());
      }
      #endif
   }

   #else
   /*
   * Compute all potential energy contributions.
   */
   void Simulation::computePotentialEnergies()
   {
      pairPotential().computeEnergy();
      bondPotential().computeEnergy();
      #ifdef INTER_ANGLE
      if (nAngleType_) {
         anglePotential().computeEnergy();
      }
      #endif
      #ifdef INTER_DIHEDRAL
      if (nDihedralType_) {
         dihedralPotential().computeEnergy();
      }
      #endif
      #ifdef INTER_EXTERNAL
      if (hasExternal_) {
         externalPotential().computeEnergy();
      }
      #endif
   }
   #endif

   /*
   * Methods computePotentialEnergies and computeStress rely on use
   * of lazy (as-needed) evaluation within the compute methods of
   * each of the potential energy classes: Each potential energy or
   * virial contribution is actually computed only if it has been
   * unset since the last evaluation. Lazy evaluation is implemented
   * explicitly in the Simulation class only for the kinetic energy
   * and kinetic stress, which are stored as members of Simulation.
   */

   /*
   * Return sum of potential energy contributions.
   */
   double Simulation::potentialEnergy() const
   {
      // Note: Pointers used here because ...Potential() accessors
      // return non-const references, which violate the method const.
      double energy = 0.0;
      energy += pairPotentialPtr_->energy();
      energy += bondPotentialPtr_->energy();
      #ifdef INTER_ANGLE
      if (nAngleType_) {
         energy += anglePotentialPtr_->energy();
      }
      #endif
      #ifdef INTER_DIHEDRAL
      if (nDihedralType_) {
         energy += dihedralPotentialPtr_->energy();
      }
      #endif
      #ifdef INTER_EXTERNAL
      if (hasExternal_) {
         energy += externalPotentialPtr_->energy();
      }
      #endif
      return energy;
   }

   /*
   * Mark all potential energies as unknown.
   */
   void Simulation::unsetPotentialEnergies()
   {
      pairPotential().unsetEnergy();
      bondPotential().unsetEnergy();
      #ifdef INTER_ANGLE
      if (nAngleType_) {
         anglePotential().unsetEnergy();
      }
      #endif
      #ifdef INTER_DIHEDRAL
      if (nDihedralType_) {
         dihedralPotential().unsetEnergy();
      }
      #endif
      #ifdef INTER_EXTERNAL
      externalPotential().unsetEnergy();
      #endif
   }

   // --- Virial Stress Methods ----------------------------------------

   #ifdef UTIL_MPI
   /*
   * Compute all virial stress contributions.
   */
   void Simulation::computeVirialStress()
   {
      pairPotential().computeStress(domain_.communicator());
      bondPotential().computeStress(domain_.communicator());
      #ifdef INTER_ANGLE
      if (nAngleType_) {
         anglePotential().computeStress(domain_.communicator());
      }
      #endif
      #ifdef INTER_DIHEDRAL
      if (nDihedralType_) {
         dihedralPotential().computeStress(domain_.communicator());
      }
      #endif
   }
   #else
   /*
   * Return stored value of virial stress (call only on master).
   */
   void Simulation::computeVirialStress()
   {
      pairPotential().computeStress();
      bondPotential().computeStress();
      #ifdef INTER_ANGLE
      if (nAngleType_) {
         anglePotential().computeStress();
      }
      #endif
      #ifdef INTER_DIHEDRAL
      if (nDihedralType_) {
         dihedralPotential().computeStress();
      }
      #endif
   }
   #endif

   /*
   * Return total virial stress.
   */
   Tensor Simulation::virialStress() const
   {
      // Note: Pointers used here because ...Potential() accessors
      // return non-const references, which violate the method const.
      Tensor stress;
      stress.zero();
      stress += pairPotentialPtr_->stress();
      stress += bondPotentialPtr_->stress();
      #ifdef INTER_ANGLE
      if (nAngleType_) {
         stress += anglePotentialPtr_->stress();
      }
      #endif
      #ifdef INTER_DIHEDRAL
      if (nDihedralType_) {
         stress += dihedralPotentialPtr_->stress();
      }
      #endif
      return stress;
   }

   /*
   * Return total virial pressure contribution.
   */
   double Simulation::virialPressure() const
   {
      double pressure;
      pressure = 0;
      pressure += pairPotentialPtr_->pressure();
      pressure += bondPotentialPtr_->pressure();
      #ifdef INTER_ANGLE
      if (nAngleType_) {
         pressure += anglePotentialPtr_->pressure();
      }
      #endif
      #ifdef INTER_DIHEDRAL
      if (nDihedralType_) {
         pressure += dihedralPotentialPtr_->pressure();
      }
      #endif
      return pressure;
   }

   #ifdef UTIL_MPI
   /*
   * Compute all pair energy contributions.
   */
   void Simulation::computePairEnergies()
   {  pairPotential().computePairEnergies(domain_.communicator()); }
   #else
   /*
   * Compute all pair energy contributions.
   */
   void Simulation::computePairEnergies()
   {  pairPotential().computePairEnergies(); }
   #endif

   /*
   * Return pair energies contributions.
   */
   DMatrix<double> Simulation::pairEnergies() const
   {
      DMatrix<double> pairEnergies;
      pairEnergies.allocate(nAtomType_, nAtomType_);
      pairEnergies = pairPotentialPtr_->pairEnergies();
      return pairEnergies;
   }

   /*
   * Mark all potential energies as unknown.
   */
   void Simulation::unsetVirialStress()
   {
      pairPotential().unsetStress();
      bondPotential().unsetStress();
      #ifdef INTER_ANGLE
      if (nAngleType_) {
         anglePotential().unsetStress();
      }
      #endif
      #ifdef INTER_DIHEDRAL
      if (nDihedralType_) {
         dihedralPotential().unsetStress();
      }
      #endif
   }

   // --- ConfigIo Accessors -------------------------------------------
   
   /**
   * Return the ConfigIo (create default if necessary).
   */
   ConfigIo& Simulation::configIo() 
   {
      if (configIoPtr_ == 0) {
         configIoPtr_ = new DdMdConfigIo(*this);
         configIoPtr_->initialize();
      }
      return *configIoPtr_;
   }

   /**
   * Return a SerialConfigIo (create if necessary).
   */
   SerializeConfigIo& Simulation::serializeConfigIo() 
   {
      if (serializeConfigIoPtr_ == 0) {
         serializeConfigIoPtr_ = new SerializeConfigIo(*this);
         serializeConfigIoPtr_->initialize();
      }
      return *serializeConfigIoPtr_;
   }

   // --- Config File Read and Write -----------------------------------

   /*
   * Read configuration file on master and distribute atoms.
   */
   void Simulation::readConfig(const std::string& filename)
   {
      std::ifstream inputFile;
      if (domain_.isMaster()) {
         fileMaster().openInputFile(filename, inputFile);
      }
      configIo().readConfig(inputFile, maskedPairPolicy_);
      exchanger_.exchange();
      if (domain_.isMaster()) {
         inputFile.close();
      }
   }

   /*
   * Write configuration file on master.
   */
   void Simulation::writeConfig(const std::string& filename)
   {
      std::ofstream outputFile;
      if (domain_.isMaster()) {
         fileMaster().openOutputFile(filename, outputFile);
      }
      configIo().writeConfig(outputFile);
      if (domain_.isMaster()) {
         outputFile.close();
      }
   }

   // --- Potential Factories and Styles -------------------------------

   #ifndef DDMD_NOPAIR
   /*
   * Return the PairFactory by reference.
   */
   Factory<PairPotential>& Simulation::pairFactory()
   {
      if (!pairFactoryPtr_) {
         pairFactoryPtr_ = new PairFactory(*this);
      }
      assert(pairFactoryPtr_);
      return *pairFactoryPtr_;
   }

   /*
   * Get the pair style string.
   */
   std::string Simulation::pairStyle() const
   {  return pairStyle_;  }
   #endif

   /*
   * Return the BondFactory by reference.
   */
   Factory<BondPotential>& Simulation::bondFactory()
   {
      if (!bondFactoryPtr_) {
         bondFactoryPtr_ = new BondFactory(*this);
      }
      assert(bondFactoryPtr_);
      return *bondFactoryPtr_;
   }

   /*
   * Get the bond style string.
   */
   std::string Simulation::bondStyle() const
   {  return bondStyle_;  }

   #ifdef INTER_ANGLE
   /*
   * Return the AngleFactory by reference.
   */
   Factory<AnglePotential>& Simulation::angleFactory()
   {
      if (angleFactoryPtr_ == 0) {
         angleFactoryPtr_ = new AngleFactory(*this);
      }
      assert(angleFactoryPtr_);
      return *angleFactoryPtr_;
   }

   /*
   * Get the angle style string.
   */
   std::string Simulation::angleStyle() const
   {  return angleStyle_;  }
   #endif

   #ifdef INTER_DIHEDRAL
   /*
   * Return the DihedralFactory by reference.
   */
   Factory<DihedralPotential>& Simulation::dihedralFactory()
   {
      if (dihedralFactoryPtr_ == 0) {
         dihedralFactoryPtr_ = new DihedralFactory(*this);
      }
      assert(dihedralFactoryPtr_);
      return *dihedralFactoryPtr_;
   }

   /*
   * Get the dihedral style string.
   */
   std::string Simulation::dihedralStyle() const
   {  return dihedralStyle_;  }
   #endif

   #ifdef INTER_EXTERNAL
   /*
   * Return the ExternalFactory by reference.
   */
   Factory<ExternalPotential>& Simulation::externalFactory()
   {
      if (externalFactoryPtr_ == 0) {
         externalFactoryPtr_ = new ExternalFactory(*this);
      }
      assert(externalFactoryPtr_);
      return *externalFactoryPtr_;
   }

   /*
   * Get the external style string.
   */
   std::string Simulation::externalStyle() const
   {  return externalStyle_;  }
   #endif

   // --- Integrator and ConfigIo Management ---------------------------

   /*
   * Return the IntegratorFactory by reference.
   */
   Factory<Integrator>& Simulation::integratorFactory()
   {
      if (integratorFactoryPtr_ == 0) {
         integratorFactoryPtr_ = new IntegratorFactory(*this);
      }
      assert(integratorFactoryPtr_);
      return *integratorFactoryPtr_;
   }

   /*
   * Return the ConfigIoFactory by reference.
   */
   Factory<ConfigIo>& Simulation::configIoFactory()
   {
      if (configIoFactoryPtr_ == 0) {
         configIoFactoryPtr_ = new ConfigIoFactory(*this);
      }
      assert(configIoFactoryPtr_);
      return *configIoFactoryPtr_;
   }

   /*
   * Set the ConfigIo, identified by subclass name.
   */
   void Simulation::setConfigIo(std::string& classname)
   {
      ConfigIo* ptr = configIoFactory().factory(classname);
      if (!ptr) {
         UTIL_THROW("Unrecognized ConfigIo subclass name");
      }
      if (configIoPtr_) {
         delete configIoPtr_;
      }
      configIoPtr_ = ptr;
      configIoPtr_->initialize();
   }

   // --- Validation ---------------------------------------------------

   /*
   * Return true if this Simulation is valid, or throw an Exception.
   */
   bool Simulation::isValid()
   {
      #ifdef UTIL_MPI
      atomStorage_.isValid(domain_.communicator());
      #else
      atomStorage_.isValid();
      #endif

      // Determine if there are any ghosts on any processor
      int nGhost = atomStorage_.nGhost();
      int nGhostAll = 0;
      #ifdef UTIL_MPI
      domain_.communicator().Reduce(&nGhost, &nGhostAll, 1,
                                    MPI::INT, MPI::SUM, 0);
      bcast(domain_.communicator(), nGhostAll, 0);
      #else
      nGhostAll = nGhost;
      #endif
      bool hasGhosts = bool(nGhostAll);

      // Test Group storage containers
      #ifdef UTIL_MPI
      bondStorage_.isValid(atomStorage_, domain_.communicator(), hasGhosts);
      #ifdef INTER_ANGLE
      if (nAngleType_) {
         angleStorage_.isValid(atomStorage_, domain_.communicator(), hasGhosts);
      }
      #endif
      #ifdef INTER_DIHEDRAL
      if (nDihedralType_) {
         dihedralStorage_.isValid(atomStorage_, domain_.communicator(),
                                  hasGhosts);
      }
      #endif
      #else // ifdef UTIL_MPI
      bondStorage_.isValid(atomStorage_, hasGhosts);
      #ifdef INTER_ANGLE
      if (nAngleType_) {
         angleStorage_.isValid(atomStorage_, hasGhosts);
      }
      #endif
      #ifdef INTER_DIHEDRAL
      if (nDihedralType_) {
         dihedralStorage_.isValid(atomStorage_, hasGhosts);
      }
      #endif
      #endif // ifdef UTIL_MPI

      // Test consistency of computed potential energies and stresses
      #ifdef UTIL_MPI
      pairPotential().isValid(domain_.communicator());
      bondPotential().isValid(domain_.communicator());
      #ifdef INTER_ANGLE
      anglePotential().isValid(domain_.communicator());
      #endif
      #ifdef INTER_DIHEDRAL
      dihedralPotential().isValid(domain_.communicator());
      #endif
      #endif // ifdef UTIL_MPI

      return true;
   }

}
#endif
