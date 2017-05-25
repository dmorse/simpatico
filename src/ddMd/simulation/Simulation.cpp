/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Simulation.h"
#include <ddMd/storage/AtomIterator.h>
#include <ddMd/storage/GhostIterator.h>
#include <ddMd/storage/GroupStorage.tpp>
#include <ddMd/integrators/Integrator.h>
#include <ddMd/integrators/IntegratorFactory.h>
#include <ddMd/configIos/ConfigIo.h>
#include <ddMd/configIos/ConfigIoFactory.h>
#include <ddMd/configIos/DdMdConfigIo.h>
#include <ddMd/configIos/SerializeConfigIo.h>
#include <ddMd/analyzers/AnalyzerManager.h>
#ifdef DDMD_MODIFIERS
#include <ddMd/modifiers/ModifierManager.h>
#endif

#include <ddMd/potentials/pair/PairPotential.h>
#include <ddMd/potentials/pair/PairFactory.h>
#ifdef SIMP_BOND
#include <ddMd/potentials/bond/BondPotential.h>
#include <ddMd/potentials/bond/BondFactory.h>
#endif
#ifdef SIMP_ANGLE
#include <ddMd/potentials/angle/AnglePotential.h>
#include <ddMd/potentials/angle/AngleFactory.h>
#endif
#ifdef SIMP_DIHEDRAL
#include <ddMd/potentials/dihedral/DihedralPotential.h>
#include <ddMd/potentials/dihedral/DihedralFactory.h>
#endif
#ifdef SIMP_EXTERNAL
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
#include <util/misc/Memory.h>
#include <util/mpi/MpiSendRecv.h>
#include <util/misc/Timer.h>
#include <util/misc/Bit.h>
#include <util/misc/initStatic.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>
#include <util/format/Str.h>
#include <util/misc/ioUtil.h>

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
      #ifdef SIMP_BOND
      bondStorage_(),
      #endif
      #ifdef SIMP_ANGLE
      angleStorage_(),
      #endif
      #ifdef SIMP_DIHEDRAL
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
      #ifdef SIMP_BOND
      bondPotentialPtr_(0),
      #endif
      #ifdef SIMP_ANGLE
      anglePotentialPtr_(0),
      #endif
      #ifdef SIMP_DIHEDRAL
      dihedralPotentialPtr_(0),
      #endif
      #ifdef SIMP_EXTERNAL
      externalPotentialPtr_(0),
      #endif
      integratorPtr_(0),
      energyEnsemblePtr_(0),
      boundaryEnsemblePtr_(0),
      fileMasterPtr_(0),
      configIoPtr_(0),
      serializeConfigIoPtr_(0),
      #ifdef DDMD_MODIFIERS
      modifierManagerPtr_(0),
      #endif
      analyzerManagerPtr_(0),
      pairFactoryPtr_(0),
      #ifdef SIMP_BOND
      bondFactoryPtr_(0),
      #endif
      #ifdef SIMP_ANGLE
      angleFactoryPtr_(0),
      #endif
      #ifdef SIMP_DIHEDRAL
      dihedralFactoryPtr_(0),
      #endif
      #ifdef SIMP_EXTERNAL
      externalFactoryPtr_(0),
      #endif
      integratorFactoryPtr_(0),
      configIoFactoryPtr_(0),
      pairStyle_(),
      #ifdef SIMP_BOND
      bondStyle_(),
      #endif
      #ifdef SIMP_ANGLE
      angleStyle_(),
      #endif
      #ifdef SIMP_DIHEDRAL
      dihedralStyle_(),
      #endif
      #ifdef SIMP_EXTERNAL
      externalStyle_(),
      #endif
      nAtomType_(0),
      #ifdef SIMP_BOND
      nBondType_(0),
      #endif
      #ifdef SIMP_ANGLE
      nAngleType_(0),
      #endif
      #ifdef SIMP_DIHEDRAL
      nDihedralType_(0),
      #endif
      #ifdef SIMP_EXTERNAL
      hasExternal_(false),
      #endif
      hasAtomContext_(false),
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
      exchanger_.associate(domain_, boundary_, atomStorage_, buffer_);
      atomStorage_.associate(domain_, boundary_, buffer_);
      #ifdef SIMP_BOND
      bondStorage_.associate(domain_, atomStorage_, buffer_);
      #endif
      #ifdef SIMP_ANGLE
      angleStorage_.associate(domain_, atomStorage_, buffer_);
      #endif
      #ifdef SIMP_DIHEDRAL
      dihedralStorage_.associate(domain_, atomStorage_, buffer_);
      #endif

      fileMasterPtr_ = new FileMaster;
      energyEnsemblePtr_ = new EnergyEnsemble;
      boundaryEnsemblePtr_ = new BoundaryEnsemble;
      #ifdef DDMD_MODIFIERS
      modifierManagerPtr_ = new ModifierManager(*this);
      #endif
      analyzerManagerPtr_ = new AnalyzerManager(*this);
   }

   /*
   * Destructor.
   */
   Simulation::~Simulation()
   {
      if (pairFactoryPtr_) {
         delete pairFactoryPtr_;
      }
      if (pairPotentialPtr_) {
         delete pairPotentialPtr_;
      }
      #ifdef SIMP_BOND
      if (bondFactoryPtr_) {
         delete bondFactoryPtr_;
      }
      if (bondPotentialPtr_) {
         delete bondPotentialPtr_;
      }
      #endif
      #ifdef SIMP_ANGLE
      if (angleFactoryPtr_) {
         delete angleFactoryPtr_;
      }
      if (anglePotentialPtr_) {
         delete anglePotentialPtr_;
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (dihedralFactoryPtr_) {
         delete dihedralFactoryPtr_;
      }
      if (dihedralPotentialPtr_) {
         delete dihedralPotentialPtr_;
      }
      #endif
      #ifdef SIMP_EXTERNAL
      if (externalFactoryPtr_) {
         delete externalFactoryPtr_;
      }
      if (externalPotentialPtr_) {
         delete externalPotentialPtr_;
      }
      #endif
      if (configIoFactoryPtr_) {
         delete configIoFactoryPtr_;
      }
      if (configIoPtr_) {
         delete configIoPtr_;
      }
      if (serializeConfigIoPtr_) {
         delete serializeConfigIoPtr_;
      }
      if (integratorFactoryPtr_) {
         delete integratorFactoryPtr_;
      }
      if (integratorPtr_) {
         delete integratorPtr_;
      }
      if (analyzerManagerPtr_) {
         delete analyzerManagerPtr_;
      }
      #ifdef DDMD_MODIFIERS
      if (modifierManagerPtr_) {
         delete modifierManagerPtr_;
      }
      #endif
      if (fileMasterPtr_) {
         delete fileMasterPtr_;
      }
      if (energyEnsemblePtr_) {
         delete energyEnsemblePtr_;
      }
      if (boundaryEnsemblePtr_) {
         delete boundaryEnsemblePtr_;
      }

      #ifdef UTIL_MPI
      if (logFile_.is_open()) {
         logFile_.close();
      }
      #endif
   }

   /*
   * Process command line options.
   */
   void Simulation::setOptions(int argc, char * const * argv)
   {
      bool sFlag = false; // split communicator
      bool eFlag = false; // echo
      bool pFlag = false; // param file name
      bool rFlag = false; // restart file name
      bool cFlag = false; // command file name
      bool iFlag = false; // input prefix
      bool oFlag = false; // output prefix
      char* sArg = 0;
      char* rArg = 0;
      char* pArg = 0;
      char* cArg = 0;
      char* iArg = 0;
      char* oArg = 0;
      int  nSystem = 1;

      // Read command-line arguments
      int c;
      opterr = 0;
      while ((c = getopt(argc, argv, "es:p:r:c:i:o:")) != -1) {
         switch (c) {
         case 'e': // echo parameters
            eFlag = true;
            break;
         case 's': // split communicator
            sFlag = true;
            sArg  = optarg;
            nSystem = atoi(sArg);
            break;
         case 'p': // parameter file
            pFlag = true;
            pArg  = optarg;
            break;
         case 'r': // restart file
            rFlag = true;
            rArg  = optarg;
            break;
         case 'c': // command file
            cFlag = true;
            cArg  = optarg;
            break;
         case 'i': // input prefix
            iFlag = true;
            iArg  = optarg;
            break;
         case 'o': // output prefix
            oFlag = true;
            oArg  = optarg;
            break;
         case '?':
            Log::file() << "Unknown option -" << optopt << std::endl;
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

      // If option -p, set parameter file name
      if (pFlag) {
         if (rFlag) {
            UTIL_THROW("Cannot have both parameter and restart files");
         }
         fileMaster().setParamFileName(std::string(pArg));
      }

      // If option -c, set command file name
      if (cFlag) {
         fileMaster().setCommandFileName(std::string(cArg));
      }

      // If option -i, set path prefix for input files
      if (iFlag) {
         fileMaster().setInputPrefix(std::string(iArg));
      }

      // If option -o, set path prefix for output files
      if (oFlag) {
         fileMaster().setOutputPrefix(std::string(oArg));
      }


      // If option -r, load state from a restart file.
      if (rFlag) {
         if (isIoProcessor()) {
            Log::file() << "Begin reading restart, base file name "
                        << std::string(rArg) << std::endl;
         }
         load(std::string(rArg));
         if (isIoProcessor()) {
            Log::file() << "Done reading restart file" << std::endl;
            Log::file() << std::endl;
         }
      }

   }

   /*
   * Read parameter block, including begin and end.
   */
   void Simulation::readParam(std::istream& in)
   {
      if (!isRestarting_) {
         readBegin(in, className().c_str());
         readParameters(in);
         readEnd(in);
      }

      /*
      * Note that if isRestarting_ is set to true, this function does 
      * nothing and returns immediately. If the ddSim program is 
      * invoked with a -r option, the setOptions() function sets 
      * isRestarting_ to true and then reads the restart file whose
      * name is specified as a command line argument. Because the 
      * main ddSim program calls setOptions() before readParam(), 
      * and because the restart file contains all the information 
      * that would otherwise be given in a parameter file, there 
      * is thus no need to read a parameter file during a restart. 
      */
   }

   /*
   *  Read parameters from default parameter file.
   */
   void Simulation::readParam()
   {
      if (!isRestarting_) {  
         readParam(fileMaster().paramFile()); 
      }
      /// See comment about restarting in readParam(std::istream&)
   }

   /*
   * Read parameters, allocate memory and initialize.
   */
   void Simulation::readParameters(std::istream& in)
   {
      if (isInitialized_) {
         UTIL_THROW("Error: Called Simulation::readParameters when isInitialized");
      }

      // Read Domain and FileMaster
      readParamComposite(in, domain_);
      readFileMaster(in);

      // Read numbers of types
      read<int>(in, "nAtomType", nAtomType_);
      #ifdef SIMP_BOND
      nBondType_ = 0;
      readOptional<int>(in, "nBondType", nBondType_); 
      if (nBondType_) {
         exchanger_.addGroupExchanger(bondStorage_);
      }
      #endif
      #ifdef SIMP_ANGLE
      nAngleType_ = 0;
      readOptional<int>(in, "nAngleType", nAngleType_); 
      if (nAngleType_) {
         exchanger_.addGroupExchanger(angleStorage_);
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      nDihedralType_ = 0;
      readOptional<int>(in, "nDihedralType", nDihedralType_); 
      if (nDihedralType_) {
         exchanger_.addGroupExchanger(dihedralStorage_);
      }
      #endif
      #ifdef SIMP_EXTERNAL
      hasExternal_ = false;
      readOptional<bool>(in, "hasExternal", hasExternal_); 
      #endif

      hasAtomContext_ = false;
      readOptional<bool>(in, "hasAtomContext", hasAtomContext_); 
      Atom::setHasAtomContext(hasAtomContext_);

      // Read array of atom type descriptors
      atomTypes_.allocate(nAtomType_);
      for (int i = 0; i < nAtomType_; ++i) {
         atomTypes_[i].setId(i);
      }
      readDArray<AtomType>(in, "atomTypes", atomTypes_, nAtomType_);

      // Read storage capacities
      readParamComposite(in, atomStorage_);
      #ifdef SIMP_BOND
      if (nBondType_) {
         readParamComposite(in, bondStorage_);
      }
      #endif
      #ifdef SIMP_ANGLE
      if (nAngleType_) {
         readParamComposite(in, angleStorage_);
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (nDihedralType_) {
         readParamComposite(in, dihedralStorage_);
      }
      #endif

      readParamComposite(in, buffer_);
      readPotentialStyles(in);

      // Create and read potential energy classes

      // Pair Potential
      assert(pairPotentialPtr_ == 0);
      pairPotentialPtr_ = pairFactory().factory(pairStyle());
      if (!pairPotentialPtr_) {
         UTIL_THROW("Unknown pairStyle");
      }
      pairPotential().setReverseUpdateFlag(reverseUpdateFlag_);
      readParamComposite(in, *pairPotentialPtr_);

      #ifdef SIMP_BOND
      // Bond Potential
      assert(bondPotentialPtr_ == 0);
      if (nBondType_) {
         bondPotentialPtr_ = bondFactory().factory(bondStyle());
         if (!bondPotentialPtr_) {
            UTIL_THROW("Unknown bondStyle");
         }
         readParamComposite(in, *bondPotentialPtr_);
      }
      #endif

      #ifdef SIMP_ANGLE
      // Angle potential
      assert(anglePotentialPtr_ == 0);
      if (nAngleType_) {
         anglePotentialPtr_ = angleFactory().factory(angleStyle());
         if (!anglePotentialPtr_) {
            UTIL_THROW("Unknown angleStyle");
         }
         readParamComposite(in, *anglePotentialPtr_);
      }
      #endif

      #ifdef SIMP_DIHEDRAL
      // Dihedral potential
      assert(dihedralPotentialPtr_ == 0);
      if (nDihedralType_) {
         dihedralPotentialPtr_ = dihedralFactory().factory(dihedralStyle());
         if (!dihedralPotentialPtr_) {
            UTIL_THROW("Unknown dihedralStyle");
         }
         readParamComposite(in, *dihedralPotentialPtr_);
      }
      #endif

      #ifdef SIMP_EXTERNAL
      // External potential
      if (hasExternal_) {
         assert(externalPotentialPtr_ == 0);
         externalPotentialPtr_ = externalFactory().factory(externalStyle());
         if (!externalPotentialPtr_) {
            UTIL_THROW("Unknown externalStyle");
         }
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
      #ifdef DDMD_MODIFIERS
      readParamCompositeOptional(in, *modifierManagerPtr_);
      #endif
      readParamComposite(in, random_);
      readParamComposite(in, *analyzerManagerPtr_);

      // Finished reading parameter file. Now finish initialization:

      exchanger_.setPairCutoff(pairPotential().cutoff());
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
      #ifdef SIMP_BOND
      if (nBondType_) {
         void (BondStorage::*memberPtr)() = &BondStorage::unsetNTotal;
         exchangeSignal().addObserver(bondStorage_, memberPtr);
      }
      #endif
      #ifdef SIMP_ANGLE
      if (nAngleType_) {
         void (AngleStorage::*memberPtr)() = &AngleStorage::unsetNTotal;
         exchangeSignal().addObserver(angleStorage_, memberPtr);
      }
      #endif
      #ifdef SIMP_DIHEDRAL
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
      #ifdef SIMP_BOND
      nBondType_ = 0;
      loadParameter<int>(ar, "nBondType", nBondType_, false); // optional
      if (nBondType_) {
         exchanger_.addGroupExchanger(bondStorage_);
      }
      #endif
      #ifdef SIMP_ANGLE
      nAngleType_ = 0;
      loadParameter<int>(ar, "nAngleType", nAngleType_, false); // opt
      if (nAngleType_) {
         exchanger_.addGroupExchanger(angleStorage_);
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      nDihedralType_ = 0;
      loadParameter<int>(ar, "nDihedralType", nDihedralType_, false); // opt
      if (nDihedralType_) {
         exchanger_.addGroupExchanger(dihedralStorage_);
      }
      #endif
      #ifdef SIMP_EXTERNAL
      hasExternal_ = false;
      loadParameter<bool>(ar, "hasExternal", hasExternal_, false); // opt
      #endif

      hasAtomContext_ = false;
      loadParameter<bool>(ar, "hasAtomContext", hasAtomContext_, false); // opt
      Atom::setHasAtomContext(hasAtomContext_);

      atomTypes_.allocate(nAtomType_);
      for (int i = 0; i < nAtomType_; ++i) {
         atomTypes_[i].setId(i);
      }
      loadDArray<AtomType>(ar, "atomTypes", atomTypes_, nAtomType_);

      // Load storage capacities
      loadParamComposite(ar, atomStorage_);
      #ifdef SIMP_BOND
      if (nBondType_) {
         loadParamComposite(ar, bondStorage_);
      }
      #endif
      #ifdef SIMP_ANGLE
      if (nAngleType_) {
         loadParamComposite(ar, angleStorage_);
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (nDihedralType_) {
         loadParamComposite(ar, dihedralStorage_);
      }
      #endif
      loadParamComposite(ar, buffer_);

      // Load potentials styles and parameters
      loadPotentialStyles(ar);

      // Pair Potential
      assert(pairPotentialPtr_ == 0);
      pairPotentialPtr_ = pairFactory().factory(pairStyle());
      if (!pairPotentialPtr_) {
         UTIL_THROW("Unknown pairStyle");
      }
      loadParamComposite(ar, *pairPotentialPtr_);
      pairPotential().setReverseUpdateFlag(reverseUpdateFlag_);

      #ifdef SIMP_BOND
      // Bond Potential
      assert(bondPotentialPtr_ == 0);
      if (nBondType_) {
         bondPotentialPtr_ = bondFactory().factory(bondStyle());
         if (!bondPotentialPtr_) {
            UTIL_THROW("Unknown bondStyle");
         }
         loadParamComposite(ar, *bondPotentialPtr_);
      }
      #endif

      #ifdef SIMP_ANGLE
      // Angle potential
      assert(anglePotentialPtr_ == 0);
      if (nAngleType_) {
         anglePotentialPtr_ = angleFactory().factory(angleStyle());
         if (!anglePotentialPtr_) {
            UTIL_THROW("Unknown angleStyle");
         }
         loadParamComposite(ar, *anglePotentialPtr_);
      }
      #endif

      #ifdef SIMP_DIHEDRAL
      // Dihedral potential
      assert(dihedralPotentialPtr_ == 0);
      if (nDihedralType_) {
         dihedralPotentialPtr_ = dihedralFactory().factory(dihedralStyle());
         if (!dihedralPotentialPtr_) {
            UTIL_THROW("Unknown dihedralStyle");
         }
         loadParamComposite(ar, *dihedralPotentialPtr_);
      }
      #endif

      #ifdef SIMP_EXTERNAL
      // External potential
      assert(externalPotentialPtr_ == 0);
      if (hasExternal_) {
         externalPotentialPtr_ = externalFactory().factory(externalStyle());
         if (!externalPotentialPtr_) {
            UTIL_THROW("Unknown externalStyle");
         }
         loadParamComposite(ar, *externalPotentialPtr_);
      }
      #endif

      loadEnsembles(ar);

      // Integrator
      assert(integratorPtr_ == 0);
      std::string className;
      integratorPtr_ =
         integratorFactory().loadObject(ar, *this, className);
      if (!integratorPtr_) {
         std::string msg("Unknown Integrator subclass name: ");
         msg += className;
         UTIL_THROW("msg.c_str()");
      }
      #ifdef DDMD_MODIFIERS
      loadParamCompositeOptional(ar, *modifierManagerPtr_);
      #endif
      loadParamComposite(ar, random_);
      loadParamComposite(ar, *analyzerManagerPtr_);

      // Finished loading data from archive. Now finish initialization:

      exchanger_.setPairCutoff(pairPotential().cutoff());
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
      #ifdef SIMP_BOND
      if (nBondType_) {
         void (BondStorage::*memberPtr)() = &BondStorage::unsetNTotal;
         exchangeSignal().addObserver(bondStorage_, memberPtr);
      }
      #endif
      #ifdef SIMP_ANGLE
      if (nAngleType_) {
         void (AngleStorage::*memberPtr)() = &AngleStorage::unsetNTotal;
         exchangeSignal().addObserver(angleStorage_, memberPtr);
      }
      #endif
      #ifdef SIMP_DIHEDRAL
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
         std::ios_base::openmode mode = std::ios_base::in | std::ios_base::binary;
         fileMaster().openRestartIFile(filename, ar.file(), mode);
      }
      // ParamComposite::load() calls Simulation::loadParameters()
      load(ar);
      if (isIoProcessor()) {
         ar.file().close();
      }

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
      #ifdef SIMP_BOND
      Parameter::saveOptional(ar, nBondType_, (bool)nBondType_);
      #endif
      #ifdef SIMP_ANGLE
      Parameter::saveOptional(ar, nAngleType_, (bool)nAngleType_);
      #endif
      #ifdef SIMP_DIHEDRAL
      Parameter::saveOptional(ar, nDihedralType_, (bool)nDihedralType_);
      #endif
      #ifdef SIMP_EXTERNAL
      Parameter::saveOptional(ar, hasExternal_, hasExternal_);
      #endif
      Parameter::saveOptional(ar, hasAtomContext_, hasAtomContext_);
      ar << atomTypes_;

      // Read storage capacities
      atomStorage_.save(ar);
      #ifdef SIMP_BOND
      if (nBondType_) {
         bondStorage_.save(ar);
      }
      #endif
      #ifdef SIMP_ANGLE
      if (nAngleType_) {
         angleStorage_.save(ar);
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (nDihedralType_) {
         dihedralStorage_.save(ar);
      }
      #endif

      buffer_.save(ar);

      // Potential energy styles and potential classes
      savePotentialStyles(ar);
      pairPotential().save(ar);
      #ifdef SIMP_BOND
      if (nBondType_) {
         bondPotential().save(ar);
      }
      #endif
      #ifdef SIMP_ANGLE
      if (nAngleType_) {
         anglePotential().save(ar);
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (nDihedralType_) {
         dihedralPotential().save(ar);
      }
      #endif
      #ifdef SIMP_EXTERNAL
      if (hasExternal_) {
         externalPotential().save(ar);
      }
      #endif

      // Save ensembles, integrator, modifiers, random, analyzers
      saveEnsembles(ar);
      std::string name = integrator().className();
      ar << name;
      integrator().save(ar);
      #ifdef DDMD_MODIFIERS
      modifierManager().saveOptional(ar);
      #endif
      random_.save(ar);
      analyzerManager().save(ar);
   }

   /*
   * Save state to file (open file, call save(), close file).
   */
   void Simulation::save(const std::string& filename)
   {
      // Update statistics (call on all processors).
      atomStorage_.computeStatistics(domain_.communicator());
      #ifdef SIMP_BOND
      if (nBondType_) {
         bondStorage_.computeStatistics(domain_.communicator());
      }
      #endif
      #ifdef SIMP_ANGLE
      if (nAngleType_) {
         angleStorage_.computeStatistics(domain_.communicator());
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (nDihedralType_) {
         dihedralStorage_.computeStatistics(domain_.communicator());
      }
      #endif
      buffer_.computeStatistics(domain_.communicator());
      if (integratorPtr_) {
         integrator().computeStatistics();
      }

      // Save parameters (only on ioProcessor)
      Serializable::OArchive ar;
      if (isIoProcessor()) {
         std::ios_base::openmode mode = std::ios_base::out | std::ios_base::binary;
         fileMaster().openRestartOFile(filename, ar.file(), mode);
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
   {
      assert(fileMasterPtr_);
      readParamComposite(in, *fileMasterPtr_);
      assert(fileMasterPtr_->hasIoCommunicator());
   }

   /*
   * Load the FileMaster.
   */
   void Simulation::loadFileMaster(Serializable::IArchive& ar)
   {
      assert(fileMasterPtr_);
      loadParamComposite(ar, *fileMasterPtr_);
      assert(fileMasterPtr_->hasIoCommunicator());
   }

   /*
   * Save the FileMaster.
   */
   void Simulation::saveFileMaster(Serializable::OArchive& ar)
   {  fileMasterPtr_->save(ar); }

   /*
   * Read potential style strings and maskedPairPolicy.
   */
   void Simulation::readPotentialStyles(std::istream &in)
   {
      read<std::string>(in, "pairStyle", pairStyle_);
      #ifdef SIMP_BOND
      if (nBondType_) {
         read<std::string>(in, "bondStyle", bondStyle_);
      }
      #endif
      #ifdef SIMP_ANGLE
      if (nAngleType_) {
         read<std::string>(in, "angleStyle", angleStyle_);
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (nDihedralType_) {
         read<std::string>(in, "dihedralStyle", dihedralStyle_);
      }
      #endif
      #ifdef SIMP_EXTERNAL
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
      loadParameter<std::string>(ar, "pairStyle", pairStyle_);
      #ifdef SIMP_BOND
      if (nBondType_) {
         loadParameter<std::string>(ar, "bondStyle", bondStyle_);
      }
      #endif
      #ifdef SIMP_ANGLE
      if (nAngleType_) {
         loadParameter<std::string>(ar, "angleStyle", angleStyle_);
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (nDihedralType_) {
         loadParameter<std::string>(ar, "dihedralStyle", dihedralStyle_);
      }
      #endif
      #ifdef SIMP_EXTERNAL
      if (hasExternal_) {
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
      ar << pairStyle_;
      #ifdef SIMP_BOND
      if (nBondType_) {
         ar << bondStyle_;
      }
      #endif
      #ifdef SIMP_ANGLE
      if (nAngleType_) {
         ar << angleStyle_;
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (nDihedralType_) {
         ar << dihedralStyle_;
      }
      #endif
      #ifdef SIMP_EXTERNAL
      if (hasExternal_) {
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
         #else // ifndef UTIL_MPI
         getNextLine(in, line);
         Log::file() << line << std::endl;
         #endif 
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
               removeDriftVelocity();
            } else
            if (command == "SIMULATE") {
               int nStep;
               inBuffer >> nStep;
               if (domain_.isMaster()) {
                  Log::file() << std::endl;
               }
               integrator().run(nStep);
            } else
            if (command == "OUTPUT_ANALYZERS") {
               analyzerManager().output();
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
               #ifdef SIMP_BOND
               if (nBondType_) {
                  bondStorage().computeStatistics(domain_.communicator());
               }
               #endif
               #ifdef SIMP_ANGLE
               if (nAngleType_) {
                  angleStorage().computeStatistics(domain_.communicator());
               }
               #endif
               #ifdef SIMP_DIHEDRAL
               if (nDihedralType_) {
                  dihedralStorage().computeStatistics(domain_.communicator());
               }
               #endif
               pairPotential().pairList()
                              .computeStatistics(domain_.communicator());
               buffer().computeStatistics(domain_.communicator());
               int maxMemory = Memory::max(domain_.communicator());
               if (domain_.isMaster()) {
                  atomStorage().outputStatistics(Log::file());
                  #ifdef SIMP_BOND
                  if (nBondType_) {
                     bondStorage().outputStatistics(Log::file());
                  }
                  #endif
                  #ifdef SIMP_ANGLE
                  if (nAngleType_) {
                     angleStorage().outputStatistics(Log::file());
                  }
                  #endif
                  #ifdef SIMP_DIHEDRAL
                  if (nDihedralType_) {
                     dihedralStorage().outputStatistics(Log::file());
                  }
                  #endif
                  buffer().outputStatistics(Log::file());
                  pairPotential().pairList().outputStatistics(Log::file());
                  Log::file() << std::endl;
                  Log::file() << "Memory: maximum allocated for arrays = "
                              << maxMemory << " bytes" << std::endl;
                  Log::file() << std::endl;
               }

               atomStorage().clearStatistics();
               #ifdef SIMP_BOND
               if (nBondType_) {
                  bondStorage().clearStatistics();
               }
               #endif
               #ifdef SIMP_ANGLE
               if (nAngleType_) {
                  angleStorage().clearStatistics();
               }
               #endif
               #ifdef SIMP_DIHEDRAL
               if (nDihedralType_) {
                  dihedralStorage().clearStatistics();
               }
               #endif
               buffer().clearStatistics();
               pairPotential().pairList().clearStatistics();

            } else
            if (command == "CLEAR_INTEGRATOR") {
               // Clear timing, memory statistics, analyzer accumulators.
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
               if (isIoProcessor()) {
                  fileMaster().openOutputFile(filename, outputFile);
                  writeParam(outputFile);
                  outputFile.close();
               }
            } else
            if (command == "SET_CONFIG_IO") {
               // Create a ConfigIo object of specified subclass.
               // This gives file format for later reads and writes.
               std::string classname;
               inBuffer >> classname;
               setConfigIo(classname);
            } else
            if (command == "SET_INPUT_PREFIX") {
               // Set the FileMaster inputPrefix, which is used to
               // construct paths to input files.
               std::string prefix;
               inBuffer >> prefix;
               fileMaster().setInputPrefix(prefix);
            } else
            if (command == "SET_OUTPUT_PREFIX") {
               // Set the FileMaster outputPrefix, which is used to
               // construct paths to output files.
               std::string prefix;
               inBuffer >> prefix;
               fileMaster().setOutputPrefix(prefix);
            } else
            if (command == "SET_PAIR") {
               // Modify one parameter of a pair interaction.
               std::string paramName;
               int typeId1, typeId2;
               double value;
               inBuffer >> paramName >> typeId1 >> typeId2 >> value;
               pairPotential().set(paramName, typeId1, typeId2, value);
            } else
            #ifdef SIMP_BOND
            if (command == "SET_BOND") {
               if (nBondType_ == 0) {
                  UTIL_THROW("SET_BOND command with nBondType = 0");
               }
               // Modify one parameter of a bond interaction.
               std::string paramName;
               int typeId;
               double value;
               inBuffer >> paramName >> typeId >> value;
               bondPotential().set(paramName, typeId, value);
            } else
            #endif
            #ifdef SIMP_ANGLE
            if (command == "SET_ANGLE") {
               if (nAngleType_ == 0) {
                  UTIL_THROW("SET_ANGLE command with nAngleType = 0");
               }
               // Modify one parameter of an angle interaction.
               std::string paramName;
               int typeId;
               double value;
               inBuffer >> paramName >> typeId >> value;
               anglePotential().set(paramName, typeId, value);
            } else
            #endif
            #ifdef SIMP_DIHEDRAL
            if (command == "SET_DIHEDRAL") {
               if (nDihedralType_ == 0) {
                  UTIL_THROW("SET_DIHEDRAL command with nDihedralType = 0");
               }
               // Modify one parameter of a dihedral interaction.
               std::string paramName;
               int typeId;
               double value;
               inBuffer >> paramName >> typeId >> value;
               dihedralPotential().set(paramName, typeId, value);
            } else
            #endif
            if (command == "SET_GROUP") {
               setGroup(inBuffer);
            } else
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
   {  
      if (fileMaster().commandFileName().empty()) {
         UTIL_THROW("Empty command file name");
      }
      readCommands(fileMaster().commandFile()); 
   }

   /*
   * Set flag to specify if reverse communication is enabled.
   */
   void Simulation::setReverseUpdateFlag(bool reverseUpdateFlag)
   {
      reverseUpdateFlag_ = reverseUpdateFlag;
      if (pairPotentialPtr_) {
         pairPotential().setReverseUpdateFlag(reverseUpdateFlag);
      }
   }

   /*
   * Choose velocities from a Boltzmann distribution.
   */
    void Simulation::setBoltzmannVelocities(double temperature)
   {
      double mass;
      double scale;
      AtomIterator atomIter;
      int i;
      atomStorage_.begin(atomIter);
      for( ; atomIter.notEnd(); ++atomIter){
         mass = atomType(atomIter->typeId()).mass();
         assert(mass > 0);
         scale = sqrt(temperature/mass);
         for (i = 0; i < Dimension; ++i) {
            atomIter->velocity()[i] = scale*random_.gaussian();
         }
      }

      // Publish notification of change in velocities
      velocitySignal().notify();
   }

   /*
   * Remove the drift velocity
   */
   Vector Simulation::removeDriftVelocity()
   {
      Vector momentum(0.0);      // atom momentum
      Vector momentumLocal(0.0); // sum of momenta on processor
      Vector momentumTotal(0.0); // total momentum of system
      double mass;               // atom mass
      double massLocal = 0.0;    // sum of masses on processor
      double massTotal = 0.0;    // total momentum of system
      int j;

      // Calculate total momentum and mass on processor
      AtomIterator atomIter;
      atomStorage_.begin(atomIter);
      for( ; atomIter.notEnd(); ++atomIter){
         mass = atomType(atomIter->typeId()).mass();
         massLocal = massLocal + mass;
         for(j = 0; j<Dimension; ++j) {
            momentum[j] = atomIter->velocity()[j];
            momentum[j] *= mass;
         }
         momentumLocal += momentum;
      }

      // Compute total momentum and mass for system, by MPI all reduce
      domain_.communicator().Allreduce(&massLocal, &massTotal, 1,
                                       MPI::DOUBLE, MPI::SUM);
      domain_.communicator().Allreduce(&momentumLocal[0],
                                       &momentumTotal[0], Dimension,
                                       MPI::DOUBLE, MPI::SUM);

      // Subtract average drift velocity
      Vector drift = momentumTotal;
      drift /= massTotal;
      atomStorage_.begin(atomIter); 
      for( ; atomIter.notEnd(); ++atomIter) {
         atomIter->velocity() -= drift;
      }

      // Publish notification of change in velocities
      velocitySignal().notify();

      return drift;
   }

   /*
   * Set forces on all atoms to zero.
   *
   * If reverseUpdateFlag() is true, zero local and ghost
   * atom forces, otherwise only local atoms.
   */
   void Simulation::zeroForces()
   {  atomStorage_.zeroForces(reverseUpdateFlag_); }

   /*
   * Compute forces for all atoms.
   */
   void Simulation::computeForces()
   {
      zeroForces();
      pairPotential().computeForces();
      #ifdef SIMP_BOND
      if (nBondType_) {
         bondPotential().computeForces();
      }
      #endif
      #ifdef SIMP_ANGLE
      if (nAngleType_) {
         anglePotential().computeForces();
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (nDihedralType_) {
         dihedralPotential().computeForces();
      }
      #endif
      #ifdef SIMP_EXTERNAL
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
      #ifdef SIMP_BOND
      if (nBondType_) {
         bondPotential().computeForcesAndStress(domain_.communicator());
      }
      #endif
      #ifdef SIMP_ANGLE
      if (nAngleType_) {
         anglePotential().computeForcesAndStress(domain_.communicator());
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (nDihedralType_) {
         dihedralPotential().computeForcesAndStress(domain_.communicator());
      }
      #endif
      #ifdef SIMP_EXTERNAL
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
      #ifdef SIMP_BOND
      if (nBondType_) {
         bondPotential().computeEnergy(domain_.communicator());
      }
      #endif
      #ifdef SIMP_ANGLE
      if (nAngleType_) {
         anglePotential().computeEnergy(domain_.communicator());
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (nDihedralType_) {
         dihedralPotential().computeEnergy(domain_.communicator());
      }
      #endif
      #ifdef SIMP_EXTERNAL
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
      #ifdef SIMP_BOND
      if (nBondType_) {
         bondPotential().computeEnergy();
      }
      #endif
      #ifdef SIMP_ANGLE
      if (nAngleType_) {
         anglePotential().computeEnergy();
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (nDihedralType_) {
         dihedralPotential().computeEnergy();
      }
      #endif
      #ifdef SIMP_EXTERNAL
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
      double energy = 0.0;
      energy += pairPotential().energy();
      #ifdef SIMP_BOND
      if (nBondType_) {
         energy += bondPotential().energy();
      }
      #endif
      #ifdef SIMP_ANGLE
      if (nAngleType_) {
         energy += anglePotential().energy();
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (nDihedralType_) {
         energy += dihedralPotential().energy();
      }
      #endif
      #ifdef SIMP_EXTERNAL
      if (hasExternal_) {
         energy += externalPotential().energy();
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
      #ifdef SIMP_BOND
      if (nBondType_) {
         bondPotential().unsetEnergy();
      }
      #endif
      #ifdef SIMP_ANGLE
      if (nAngleType_) {
         anglePotential().unsetEnergy();
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (nDihedralType_) {
         dihedralPotential().unsetEnergy();
      }
      #endif
      #ifdef SIMP_EXTERNAL
      if (hasExternal_) {
         externalPotential().unsetEnergy();
      }
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
      #ifdef SIMP_BOND
      if (nBondType_) {
         bondPotential().computeStress(domain_.communicator());
      }
      #endif
      #ifdef SIMP_ANGLE
      if (nAngleType_) {
         anglePotential().computeStress(domain_.communicator());
      }
      #endif
      #ifdef SIMP_DIHEDRAL
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
      #ifdef SIMP_BOND
      if (nBondType_) {
         bondPotential().computeStress();
      }
      #endif
      #ifdef SIMP_ANGLE
      if (nAngleType_) {
         anglePotential().computeStress();
      }
      #endif
      #ifdef SIMP_DIHEDRAL
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
      Tensor stress;
      stress.zero();
      stress += pairPotential().stress();
      #ifdef SIMP_BOND
      if (nBondType_) {
         stress += bondPotential().stress();
      }
      #endif
      #ifdef SIMP_ANGLE
      if (nAngleType_) {
         stress += anglePotential().stress();
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (nDihedralType_) {
         stress += dihedralPotential().stress();
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
      pressure += pairPotential().pressure();
      #ifdef SIMP_BOND
      if (nBondType_) {
         pressure += bondPotential().pressure();
      }
      #endif
      #ifdef SIMP_ANGLE
      if (nAngleType_) {
         pressure += anglePotential().pressure();
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (nDihedralType_) {
         pressure += dihedralPotential().pressure();
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
      pairEnergies = pairPotential().pairEnergies();
      return pairEnergies;
   }

   /*
   * Mark all potential energies as unknown.
   */
   void Simulation::unsetVirialStress()
   {
      pairPotential().unsetStress();
      #ifdef SIMP_BOND
      if (nBondType_) {
         bondPotential().unsetStress();
      }
      #endif
      #ifdef SIMP_ANGLE
      if (nAngleType_) {
         anglePotential().unsetStress();
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (nDihedralType_) {
         dihedralPotential().unsetStress();
      }
      #endif
   }

   // --- ConfigIo Accessors -------------------------------------------

   /*
   * Return the ConfigIo (create default if necessary).
   */
   ConfigIo& Simulation::configIo()
   {
      if (configIoPtr_ == 0) {
         if (Atom::hasAtomContext()) {
            // DdMdConfigIo format with hasMolecules = true
            configIoPtr_ = new DdMdConfigIo(*this, true);
         } else {
            // DdMdConfigIo format with hasMolecules = false
            configIoPtr_ = new DdMdConfigIo(*this, false);
         }
      }
      return *configIoPtr_;
   }

   /*
   * Return a SerialConfigIo (create if necessary).
   */
   SerializeConfigIo& Simulation::serializeConfigIo()
   {
      if (serializeConfigIoPtr_ == 0) {
         serializeConfigIoPtr_ = new SerializeConfigIo(*this);
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

   #ifdef SIMP_BOND
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
   #endif

   #ifdef SIMP_ANGLE
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

   #ifdef SIMP_DIHEDRAL
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

   #ifdef SIMP_EXTERNAL
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
   }

   // --- Group Management ---------------------------------------------
   
   /*
   * Return true if this Simulation is valid, or throw an Exception.
   */
   void Simulation::setGroup(std::stringstream& inBuffer)
   {

      // Read groupId and create a corresponding bit mask
      unsigned int groupId;
      inBuffer >> groupId;
      Bit bit(groupId);

      // Read the group style
      std::string groupStyle;
      inBuffer >> groupStyle;

      AtomIterator iter;
      if (groupStyle == "delete") {
         for (atomStorage_.begin(iter) ; iter.notEnd(); ++iter) {
            bit.clear(iter->groups());
         }
      } else 
      if (groupStyle == "all") {
         for (atomStorage_.begin(iter) ; iter.notEnd(); ++iter) {
            bit.set(iter->groups());
         }
      } else 
      if (groupStyle == "atomType") {
         int typeId;
         inBuffer >> typeId;
         for (atomStorage_.begin(iter) ; iter.notEnd(); ++iter) {
            if (iter->typeId() == typeId) {
               bit.set(iter->groups());
            }
         }
      } else
      if (groupStyle == "species") {
         if (!Atom::hasAtomContext()) {
            UTIL_THROW("AtomContext is not enabled");
         }
         int speciesId;
         inBuffer >> speciesId;
         for (atomStorage_.begin(iter) ; iter.notEnd(); ++iter) {
            if (iter->context().speciesId == speciesId) {
               bit.set(iter->groups());
            }
         }
      } else {
         std::string msg = "SET_GROUP command with unknown group style ";
         msg += groupStyle;
         UTIL_THROW(msg.c_str());
      }
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
      #ifdef SIMP_BOND
      if (nBondType_) {
         bondStorage_.isValid(atomStorage_, domain_.communicator(), hasGhosts);
      }
      #endif
      #ifdef SIMP_ANGLE
      if (nAngleType_) {
         angleStorage_.isValid(atomStorage_, domain_.communicator(), hasGhosts);
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (nDihedralType_) {
         dihedralStorage_.isValid(atomStorage_, domain_.communicator(),
                                  hasGhosts);
      }
      #endif
      #else // ifdef UTIL_MPI
      #ifdef SIMP_BOND
      if (nBondType_) {
         bondStorage_.isValid(atomStorage_, hasGhosts);
      }
      #endif
      #ifdef SIMP_ANGLE
      if (nAngleType_) {
         angleStorage_.isValid(atomStorage_, hasGhosts);
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (nDihedralType_) {
         dihedralStorage_.isValid(atomStorage_, hasGhosts);
      }
      #endif
      #endif // ifdef UTIL_MPI

      // Test consistency of computed potential energies and stresses
      #ifdef UTIL_MPI
      pairPotential().isValid(domain_.communicator());
      #ifdef SIMP_BOND
      if (nBondType_) {
         bondPotential().isValid(domain_.communicator());
      }
      #endif
      #ifdef SIMP_ANGLE
      if (nAngleType_) {
         anglePotential().isValid(domain_.communicator());
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (nDihedralType_) {
         dihedralPotential().isValid(domain_.communicator());
      }
      #endif
      #ifdef SIMP_EXTERNAL
      if (hasExternal_) {
         externalPotential().isValid(domain_.communicator());
      }
      #endif
      #endif // ifdef UTIL_MPI

      return true;
   }

}
