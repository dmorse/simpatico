#ifndef DDMD_SYSTEM_CPP
#define DDMD_SYSTEM_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "System.h"
#include <ddMd/storage/AtomIterator.h>
#include <ddMd/potentials/pair/PairPotential.h>
#include <ddMd/potentials/pair/PairPotentialImpl.h>
#include <ddMd/potentials/bond/BondPotential.h>
#include <ddMd/potentials/bond/BondPotentialImpl.h>
#include <ddMd/integrator/Integrator.h>
#include <ddMd/integrator/IntegratorFactory.h>
#include <ddMd/integrator/NveIntegrator.h>
#include <ddMd/configIo/ConfigIo.h>
#include <ddMd/util/FileMaster.h>
#include <ddMd/diagnostics/DiagnosticManager.h>

#ifndef DDMD_NOPAIR
#include <ddMd/potentials/pair/PairFactory.h>
#endif
#include <ddMd/potentials/bond/BondFactory.h>
#ifdef DDMD_ANGLE
#include <ddMd/potentials/angle/AngleFactory.h>
#endif
#ifdef DDMD_DIHEDRAL
#include <ddMd/potentials/dihedral/DihedralFactory.h>
#endif
#ifdef DDMD_EXTERNAL
#include <ddMd/potentials/external/ExternalFactory.h>
#endif
#if 0
#include <ddMd/integrator/IntegratorFactory.h>
#endif
#if 0
#include <ddMd/configIos/ConfigIoFactory.h>
#endif

// namespace McMd
#include <mcMd/mcSimulation/McSimulation.h>

// namespace Util
#include <util/ensembles/EnergyEnsemble.h>
#include <util/ensembles/BoundaryEnsemble.h>
#include <util/space/Vector.h>
#include <util/space/IntVector.h>
#include <util/param/Factory.h>
#include <util/util/Log.h>
#include <util/mpi/MpiSendRecv.h>
#include <util/util/Timer.h>

#include <util/format/Int.h>
#include <util/format/Dbl.h>
#include <util/format/Str.h>
#include <util/util/ioUtil.h>
#include <util/util/initStatic.h>

#include <fstream>

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   #ifdef UTIL_MPI
   System::System(MPI::Intracomm& communicator)
   #else
   System::System() 
   #endif
    : atomStorage_(),
      bondStorage_(),
      boundary_(),
      atomTypes_(),
      domain_(),
      buffer_(),
      exchanger_(),
      random_(),
      maxBoundary_(),
      kineticEnergy_(0.0),
      #ifdef UTIL_MPI
      //communicatorPtr_(&communicator),
      #endif
      pairPotentialPtr_(0),
      bondPotentialPtr_(0),
      integratorPtr_(0),
      configIoPtr_(0),
      energyEnsemblePtr_(0),
      boundaryEnsemblePtr_(0),
      fileMasterPtr_(0),
      #ifndef DDMD_NOPAIR
      pairFactoryPtr_(0),
      #endif
      bondFactoryPtr_(0),
      #ifdef DDMD_ANGLE
      angleFactoryPtr_(0),
      #endif
      #ifdef DDMD_DIHEDRAL
      dihedralFactoryPtr_(0),
      #endif
      integratorFactoryPtr_(0),
      #if 0
      configIoFactoryPtr_(0),
      #endif
      #ifndef DDMD_NOPAIR
      pairStyle_(),
      #endif
      bondStyle_(),
      #ifdef DDMD_ANGLE
      angleStyle_(),
      #endif
      #ifdef DDMD_DIHEDRAL
      dihedralStyle_(),
      #endif
      nAtomType_(0),
      nBondType_(0),
      maskedPairPolicy_(MaskBonded)
   {
      Util::initStatic();

      #ifdef UTIL_MPI
      if (!MPI::Is_initialized()) {
         UTIL_THROW("MPI is not initialized");
      }
      Util::Vector::commitMpiType();
      Util::IntVector::commitMpiType();
      AtomType::initStatic();

      setParamCommunicator(communicator);
      domain_.setGridCommunicator(communicator);
      #endif

      // Set connections between member objects
      domain_.setBoundary(boundary_);
      exchanger_.associate(domain_, boundary_,
                           atomStorage_, bondStorage_, buffer_);

      energyEnsemblePtr_  = new EnergyEnsemble;
      boundaryEnsemblePtr_ = new BoundaryEnsemble;
     
      diagnosticManagerPtr_ = new DiagnosticManager(*this);
   }

   /*
   * Destructor.
   */
   System::~System()
   {
      if (pairPotentialPtr_) {
         delete pairPotentialPtr_;
      }
      if (bondPotentialPtr_) {
         delete bondPotentialPtr_;
      }
      if (energyEnsemblePtr_) {
         delete energyEnsemblePtr_;
      }
      if (boundaryEnsemblePtr_) {
         delete boundaryEnsemblePtr_;
      }
      if (integratorPtr_) {
         delete integratorPtr_;
      }
      if (configIoPtr_) {
         delete configIoPtr_;
      }
      if (fileMasterPtr_) {
         delete fileMasterPtr_;
      }
      #ifndef DDMD_NOPAIR
      if (pairFactoryPtr_) {
         delete pairFactoryPtr_;
      }
      #endif
      if (bondFactoryPtr_) {
         delete bondFactoryPtr_;
      }
      #ifdef DDMD_ANGLE
      if (angleFactoryPtr_) {
         delete angleFactoryPtr_;
      }
      #endif
      #ifdef DDMD_DIHEDRAL
      if (dihedralFactoryPtr_) {
         delete dihedralFactoryPtr_;
      }
      #endif
      #ifdef DDMD_EXTERNAL
      if (externalFactoryPtr_) {
         delete externalFactoryPtr_;
      }
      #endif
      if (integratorFactoryPtr_) {
         delete integratorFactoryPtr_;
      }
      #if 0
      if (configIoFactoryPtr_) {
         delete configIoFactoryPtr_;
      }
      #endif
      if (diagnosticManagerPtr_) {
         delete diagnosticManagerPtr_;
      }

      #ifdef UTIL_MPI
      //if (logFile_.is_open()) logFile_.close();
      #endif
   }

   /**
   * Read parameters, allocate memory and initialize.
   */
   void System::readParam(std::istream& in)
   {

      // Preconditions
      assert(pairPotentialPtr_ == 0);
      assert(bondPotentialPtr_ == 0);
      assert(integratorPtr_ == 0);
      assert(configIoPtr_ == 0);

      readBegin(in, "System");
      readParamComposite(in, domain_);
      readParamComposite(in, atomStorage_);
      readParamComposite(in, bondStorage_);
      readParamComposite(in, buffer_);

      // Read types
      readFileMaster(in);
      read<int>(in, "nAtomType", nAtomType_);
      read<int>(in, "nBondType", nBondType_);
      atomTypes_.allocate(nAtomType_);
      for (int i = 0; i < nAtomType_; ++i) {
         atomTypes_[i].setId(i);
      }
      readDArray<AtomType>(in, "atomTypes", atomTypes_, nAtomType_);

      readPotentialStyles(in);

      #ifndef DDMD_NOPAIR
      // Pair Potential
      pairPotentialPtr_ = pairFactory().factory(pairStyle());
      pairPotentialPtr_->setNAtomType(nAtomType_);
      readParamComposite(in, *pairPotentialPtr_);
      #endif

      // Bond Potential
      bondPotentialPtr_ = bondFactory().factory(bondStyle());
      bondPotentialPtr_->setNBondType(nBondType_);
      readParamComposite(in, *bondPotentialPtr_);

      readEnsembles(in);

      // Integrator
      std::string className;
      bool        isEnd;
      integratorPtr_ = 
         integratorFactory().readObject(in, *this, className, isEnd);
      if (!integratorPtr_) {
         std::string msg("Unknown Integrator subclass name: ");
         msg += className;
         UTIL_THROW("msg.c_str()");
      }

      readParamComposite(in, random_);
      readParamComposite(in, *diagnosticManagerPtr_);

      configIoPtr_ = new ConfigIo();             // Todo: Add factory
      configIoPtr_->associate(domain_, boundary_,
                              atomStorage_, bondStorage_, buffer_);
      readParamComposite(in, *configIoPtr_);

      exchanger_.setPairCutoff(pairPotentialPtr_->cutoff());
      exchanger_.allocate();

      readEnd(in);
   }

   /**
   * If no FileMaster exists, create and initialize one. 
   */
   void System::readFileMaster(std::istream &in)
   {
      // Create FileMaster if necessary
      if (!fileMasterPtr_) {
         fileMasterPtr_ = new FileMaster();
         readParamComposite(in, *fileMasterPtr_);
      }
   }

   /**
   * Read potential style strings and maskedPairPolicy.
   */
   void System::readPotentialStyles(std::istream &in)
   {
      #ifndef DDMD_NOPAIR
      read<std::string>(in, "pairStyle", pairStyle_);
      #endif

      if (nBondType() > 0) {
         read<std::string>(in, "bondStyle", bondStyle_);
      }

      #ifdef DDMD_ANGLE
      if (nAngleType() > 0) {
         read<std::string>(in, "angleStyle", angleStyle_);
      }
      #endif

      #ifdef DDMD_DIHEDRAL
      if (nDihedralType() > 0) {
         read<std::string>(in, "dihedralStyle", dihedralStyle_);
      }
      #endif
      #ifdef DDMD_EXTERNAL
      if (simulation().hasExternal()) {
         read<std::string>(in, "externalStyle", externalStyle_);
      }
      #endif

      read<MaskPolicy>(in, "maskedPairPolicy", maskedPairPolicy_);
   }

   /*
   * Read EnergyEnsemble and BoundaryEnsemble
   */
   void System::readEnsembles(std::istream &in)
   {
      readParamComposite(in, *energyEnsemblePtr_);
      readParamComposite(in, *boundaryEnsemblePtr_);
   }

   /*
   * Read and implement commands in an input script.
   */
   void System::readCommands(std::istream &in)
   {
      std::string   command;
      std::string   filename;
      std::ifstream inputFile;
      std::ofstream outputFile;

      #ifndef UTIL_MPI
      std::istream&     inBuffer = in;
      #else
      std::stringstream inBuffer;
      std::string       line;
      #endif

      bool readNext = true;
      while (readNext) {

         #ifdef UTIL_MPI
         if (!hasParamCommunicator() || isParamIoProcessor()) {
            getNextLine(in, line);
            Log::file() << line << std::endl;
         } 
         if (hasParamCommunicator()) {
            bcast<std::string>(domain_.communicator(), line, 0);
         }
         inBuffer.clear();
         for (unsigned i=0; i < line.size(); ++i) {
            inBuffer.put(line[i]);
         }
         #endif

         inBuffer >> command;
         //Log::file().setf(std::ios_base::left);
         //Log::file().width(15);
         //Log::file() << command;

         if (command == "READ_CONFIG") {
            inBuffer >> filename;
            readConfig(filename);
            #if 0
            if (domain_.isMaster()) {
               fileMaster().openInputFile(filename, inputFile);
            }
            configIoPtr_->readConfig(inputFile, maskedPairPolicy_);
            exchanger_.exchange();
            if (domain_.isMaster()) {
               inputFile.close();
            }
            #endif
         } else
         if (command == "THERMALIZE") {
            double temperature;
            inBuffer >> temperature;
            //Log::file() << Dbl(temperature, 15, 6) << std::endl;
            setBoltzmannVelocities(temperature);
            //removeDriftVelocity();
         } else
         if (command == "SIMULATE") {
            int endStep;
            inBuffer >> endStep;
            simulate(endStep);
         } else
         if (command == "WRITE_CONFIG") {
            inBuffer >> filename;
            writeConfig(filename);
            #if 0
            if (domain_.isMaster()) {
               fileMaster().openOutputFile(filename, outputFile);
            }
            configIoPtr_->writeConfig(outputFile);
            if (domain_.isMaster()) {
               outputFile.close();
            }
            #endif
         } else
         if (command == "WRITE_PARAM") {
            inBuffer >> filename;
            fileMaster().openOutputFile(filename, outputFile);
            //outputFile.open(filename.c_str());
            writeParam(outputFile);
            outputFile.close();
         } else
         #if 0
         if (command == "SET_CONFIG_IO") {
            std::string classname;
            inBuffer >> classname;
            Log::file() << Str(classname, 15) << std::endl;
            setConfigIo(classname);
         } else
         #endif
         if (command == "FINISH") {
            readNext = false;
         } else {
            Log::file() << "Error: Unknown command  " << std::endl;
            readNext = false;
         }

      }
   }

   /*
   * Read and implement commands from the default command file.
   */
   void System::readCommands()
   {  readCommands(fileMaster().commandFile()); }

   /*
   * Choose velocities from a Boltzmann distribution.
   */
   void System::setBoltzmannVelocities(double temperature)
   {
      double scale = sqrt(temperature);
      AtomIterator atomIter;
      int i;

      atomStorage_.begin(atomIter); 
      for( ; !atomIter.atEnd(); ++atomIter){
         for (i = 0; i < Dimension; ++i) {
            atomIter->velocity()[i] = scale*random_.gaussian();
         }
      }
   }

   /*
   * Set forces on all local atoms to zero.
   */
   void System::zeroForces()
   {
      AtomIterator atomIter;
      atomStorage_.begin(atomIter); 
      for( ; !atomIter.atEnd(); ++atomIter){
         atomIter->force().zero();
      }
   }

   /*
   * Set forces on all local atoms to zero.
   */
   void System::computeForces()
   {
      zeroForces();
      pairPotential().addForces();
      bondPotential().addForces();
   }

   /*
   * Integrate.
   */
   void System::simulate(int nStep)
   {
      // Preconditions
      assert(integratorPtr_);

      Timer timer;
      bool isMaster = bool(domain_.isMaster());
      if (isMaster) {
         Log::file() << std::endl;
      }

      integratorPtr_->setup();
      diagnosticManager().setup();

      // Main MD loop
      timer.start();
      for (int iStep_ = 0; iStep_ < nStep; ++iStep_) {

         if (Diagnostic::baseInterval > 0) {
            if (iStep_ % Diagnostic::baseInterval == 0) {
               diagnosticManager().sample(iStep_);
            }
         }

         integratorPtr_->step();
      }
      timer.stop();

      atomStorage_.computeNAtomTotal(domain_.communicator());

      if (isMaster) {

         int nAtomTot = atomStorage_.nAtomTotal();
         int nProc = 1;
         #ifdef UTIL_MPI
         nProc = domain_.communicator().Get_size();
         #endif

         // Output total time for the run
         Log::file() << std::endl;
         Log::file() << "Time Statistics" << std::endl;
         Log::file() << "nStep                " << nStep << std::endl;
         Log::file() << "run time             " << timer.time() << " sec" << std::endl;
         Log::file() << "time / nStep         " << timer.time()/double(nStep) 
   	             << " sec" << std::endl;
         Log::file() << "time / (nStep*nAtom) " 
                     << timer.time()*double(nProc)/double(nStep*nAtomTot)
                     << " sec" << std::endl;
         Log::file() << std::endl;
         Log::file() << std::endl;
  
         #if 0 
         Log::file() << "PairList Statistics" << std::endl;
         Log::file() << "maxNPair           " 
                     << pairPotential().pairList().maxNPair()
                     << std::endl;
         Log::file() << "maxNAtom           " 
                     << pairPotential().pairList().maxNAtom()
                     << std::endl;
         Log::file() << "buildCounter       " 
                     << pairPotential().pairList().buildCounter()
                     << std::endl;
         Log::file() << "steps / build      "
                     << double(nStep)/double(pairPotential().pairList().buildCounter())
                     << std::endl;
         Log::file() << std::endl;
         #endif

      }

   }

   /*
   * Calculate total kinetic energy
   *
   * Call on all processors. 
   */
   void System::computeKineticEnergy()
   {
      double energy = 0.0;
      double mass;
      int typeId;

      // Add kinetic energies of local atoms on this processor
      AtomIterator atomIter;
      atomStorage_.begin(atomIter); 
      for( ; atomIter.notEnd(); ++atomIter){
         typeId = atomIter->typeId();
         mass   = atomTypes_[typeId].mass();
         energy += mass*(atomIter->velocity().square());
      }
      energy = 0.5*energy;

      #ifdef UTIL_MPI
      // Sum values from all processors.
      domain_.communicator().Reduce(&energy, &kineticEnergy_, 1, 
                                    MPI::DOUBLE, MPI::SUM, 0);
      #else
      kineticEnergy_ = energy;
      #endif

   }

   /*
   * Return total kinetic energy (on master processor).
   * 
   * Called only on master processor.
   */
   double System::kineticEnergy()
   {  return kineticEnergy_; }

   /*
   * Compute all potential energy contributions.
   */
   void System::computePotentialEnergies() 
   {
      pairPotential().computeEnergy(domain_.communicator());
      bondPotential().computeEnergy(domain_.communicator());
   }

   /*
   * Compute all potential energy contributions.
   */
   double System::potentialEnergy() 
   {
      double energy = 0.0;
      energy += pairPotential().energy();
      energy += bondPotential().energy();
      return energy;
   }

   /*
   * Read configuration file on master and distribute atoms.
   *
   * \param filename name of configuration file.
   */
   void System::readConfig(const std::string& filename)
   {
      assert(configIoPtr_);
      std::ifstream inputFile;
      if (domain_.isMaster()) {
         fileMaster().openInputFile(filename, inputFile);
      }
      configIoPtr_->readConfig(inputFile, maskedPairPolicy_);
      exchanger_.exchange();
      if (domain_.isMaster()) {
         inputFile.close();
      }
   }

   /*
   * Write configuration file on master.
   *
   * \param filename name of configuration file.
   */
   void System::writeConfig(const std::string& filename)
   {
      assert(configIoPtr_);
      std::ofstream outputFile;
      if (domain_.isMaster()) {
         fileMaster().openOutputFile(filename, outputFile);
      }
      configIoPtr_->writeConfig(outputFile);
      if (domain_.isMaster()) {
         outputFile.close();
      }
   }

   /**
   * Determine whether an atom exchange and reneighboring is needed.
   */
   bool System::needExchange() 
   {
      // Calculate maximum square displacment among along nodes
      double maxSqDisp = atomStorage_.maxSqDisplacement();
      double maxSqDispAll;
      #if UTIL_MPI
      domain_.communicator().Reduce(&maxSqDisp, &maxSqDispAll, 1, 
                                    MPI::DOUBLE, MPI::MAX, 0);
      #else
      maxSqDispAll = maxSqDisp;
      #endif

      // Decide if maximum exceeds threshhold (on master)
      int needed = 0;
      if (domain_.isMaster()) {
         if (sqrt(maxSqDispAll) > 0.5*pairPotentialPtr_->skin()) {
            needed = 1; 
         }
      }

      #if UTIL_MPI
      // Broadcast decision to all nodes
      domain_.communicator().Bcast(&needed, 1, MPI::INT, 0);
      #endif

      return bool(needed);
   }

   #ifndef DDMD_NOPAIR
   /*
   * Return the PairFactory by reference.
   */
   Factory<PairPotential>& System::pairFactory()
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
   std::string System::pairStyle() const
   {  return pairStyle_;  }
   #endif

   /*
   * Return the BondFactory by reference.
   */
   Factory<BondPotential>& System::bondFactory()
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
   std::string System::bondStyle() const
   {  return bondStyle_;  }

   #ifdef DDMD_ANGLE
   /*
   * Return the AngleFactory by reference.
   */
   Factory<AnglePotential>& System::angleFactory()
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
   std::string System::angleStyle() const
   {  return angleStyle_;  }
   #endif

   #ifdef DDMD_DIHEDRAL
   /*
   * Return the DihedralFactory by reference.
   */
   Factory<DihedralPotential>& System::dihedralFactory()
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
   std::string System::dihedralStyle() const
   {  return dihedralStyle_;  }
   #endif

   #ifdef DDMD_EXTERNAL
   /*
   * Return the ExternalFactory by reference.
   */
   Factory<ExternalPotential>& System::externalFactory()
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
   std::string System::externalStyle() const
   {  return externalStyle_;  }
   #endif

   /*
   * Return the IntegratorFactory by reference.
   */
   Factory<Integrator>& System::integratorFactory()
   {
      if (integratorFactoryPtr_ == 0) {
         integratorFactoryPtr_ = new IntegratorFactory(*this);
      }
      assert(integratorFactoryPtr_);
      return *integratorFactoryPtr_;
   }

   #if 0
   // ConfigIoIo Management

   /*
   * Get the ConfigIo factory by reference.
   */
   Factory<ConfigIo>& System::configIoFactory()
   {
      if (!configIoFactoryPtr_) {
         configIoFactoryPtr_ = newDefaultConfigIoFactory();
      }
      return *configIoFactoryPtr_;
   }

   /*
   * Return a pointer to a new ConfigIoFactory.
   */
   Factory<ConfigIo>* System::newDefaultConfigIoFactory()
   {  return new ConfigIoFactory(*this); }

   /*
   * Set the ConfigIo, identified by subclass name.
   */
   void System::setConfigIo(std::string& classname)
   {
      if (!configIoFactoryPtr_) {
         configIoFactoryPtr_ = newDefaultConfigIoFactory();
      }
      ConfigIo* ptr = configIoFactoryPtr_->factory(classname);
      if (!ptr) {
         UTIL_THROW("Unrecognized ConfigIo subclass name");
      } 
      if (configIoPtr_) {
         delete configIoPtr_;
      }
      configIoPtr_ = ptr;
   }
   #endif

   #if 0
   /*
   * Return a pointer to a new default ConfigIo.
   */
   ConfigIo* System::newDefaultConfigIo()
   {  return new McConfigIo(*this); }
   #endif

   /**
   * Return true if this System is valid, or throw an Exception.
   */
   bool System::isValid()
   {
      atomStorage_.isValid();

      // Determine if there are any ghosts on any processor
      int nGhost = atomStorage_.nGhost();
      int nGhostAll = 0;
      #ifdef UTIL_MPI
      domain_.communicator().Reduce(&nGhost, &nGhostAll, 1, 
                                    MPI::INT, MPI::SUM, 0);
      #else
      nGhostAll = nGhost;
      #endif
      #if UTIL_MPI
      bcast(domain_.communicator(), nGhostAll, 0);
      #endif
      bool hasGhosts = bool(nGhostAll);

      bondStorage_.isValid(atomStorage_, domain_.communicator(), hasGhosts);
      return true; 

   }

}
#endif
