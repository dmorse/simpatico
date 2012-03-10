#ifndef SYSTEM_CPP
#define SYSTEM_CPP

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
#include <ddMd/integrator/NveIntegrator.h>
#include <ddMd/configIo/ConfigIo.h>

#include <ddMd/ensembles/EnergyEnsemble.h>
#include <ddMd/ensembles/BoundaryEnsemble.h>

#if 0
#include <ddMd/util/FileMaster.h>
#endif

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
#include <ddMd/configIos/ConfigIoFactory.h>
#endif

// namespace Util
#include <util/space/Vector.h>
#include <util/space/IntVector.h>
#include <util/param/Factory.h>
#include <util/util/Log.h>
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
   * Default constructor.
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
      #ifdef UTIL_MPI
      communicatorPtr_(&communicator),
      #endif
      pairPotentialPtr_(0),
      bondPotentialPtr_(0),
      integratorPtr_(0),
      configIoPtr_(0),
      energyEnsemblePtr_(0),
      boundaryEnsemblePtr_(0),
      #if 0
      fileMasterPtr_(0),
      #endif
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
      maskedPairPolicy_(MaskBonded),
      isMaster_(false)
   {
      Util::initStatic();

      #ifdef UTIL_MPI
      if (!MPI::Is_initialized()) {
         MPI::Init();
         Util::Vector::commitMpiType();
         Util::IntVector::commitMpiType();
         //Util::Pair<int>::commitMpiType();
      }

      communicatorPtr_ = &communicator;
      domain_.setGridCommunicator(communicator);
      setParamCommunicator(communicator);

      // Set directory Id in FileMaster to MPI processor rank.
      // fileMaster_.setDirectoryId(communicatorPtr_->Get_rank());

      // Set log file for processor n to a new file named "n/log"
      // Relies on initialization of FileMaster outputPrefix to "" (empty).
      // fileMaster_.openOutputFile("log", logFile_);
      // Log::setFile(logFile_);
      #endif

      // Set connections between member objects
      domain_.setBoundary(boundary_);
      exchanger_.associate(domain_, boundary_,
                           atomStorage_, bondStorage_, buffer_);

      energyEnsemblePtr_  = new EnergyEnsemble;
      boundaryEnsemblePtr_ = new BoundaryEnsemble;
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
      #if 0
      if (fileMasterPtr_) {
         delete fileMasterPtr_;
      }
      #endif
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
      #if 0
      if (configIoFactoryPtr_) {
         delete configIoFactoryPtr_;
      }
      #endif

      #ifdef UTIL_MPI
      //if (logFile_.is_open()) logFile_.close();
      MPI::Finalize();
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
      integratorPtr_ = new NveIntegrator(*this); // Todo: Add factory
      readParamComposite(in, *integratorPtr_);

      readParamComposite(in, random_);

      configIoPtr_ = new ConfigIo();             // Todo: Add factory
      configIoPtr_->associate(domain_, boundary_,
                              atomStorage_, bondStorage_, buffer_);
      readParamComposite(in, *configIoPtr_);

      exchanger_.allocate();
      exchanger_.setPairCutoff(pairPotentialPtr_->cutoff());

      readEnd(in);
   }

   #if 0
   /**
   * If no FileMaster exists, create and initialize one. 
   */
   void System::readFileMaster(std::istream &in)
   {
      // Create FileMaster if necessary
      if (!fileMasterPtr_) {
         fileMasterPtr_ = new FileMaster();
         createdFileMaster_ = true;
         readParamComposite(in, *fileMasterPtr_);
      }
   }
   #endif

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
            //Log::file() << Str(filename, 15) << std::endl;
            //fileMaster().openInputFile(filename, inputFile);
            //readConfig(inputFile);
            configIoPtr_->readConfig(filename.c_str(), maskedPairPolicy_);
            exchanger_.exchange();
            //inputFile.close();
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
            //Log::file() << Int(endStep, 15) << std::endl;
            integrate(endStep);
         } else
         if (command == "WRITE_CONFIG") {
            inBuffer >> filename;
            //Log::file() << Str(filename, 15) << std::endl;
            //fileMaster().openOutputFile(filename, outputFile);
            //writeConfig(outputFile);
            configIoPtr_->writeConfig(filename);
            //outputFile.close();
         } else
         if (command == "WRITE_PARAM") {
            inBuffer >> filename;
            //Log::file() << Str(filename, 15) << std::endl;
            //fileMaster().openOutputFile(filename, outputFile);
            outputFile.open(filename.c_str());
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
            //Log::file() << std::endl;
            readNext = false;
         } else {
            Log::file() << "Error: Unknown command  " << std::endl;
            readNext = false;
         }

      }
   }

   #if 0
   /*
   * Read and implement commands from the default command file.
   */
   void System::readCommands()
   {  readCommands(fileMaster().commandFile()); }
   #endif

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
   void System::integrate(int nStep)
   {
      // Preconditions
      assert(integratorPtr_);

      Timer timer;
      bool isMaster = bool(domain_.isMaster());
      if (isMaster) {
         Log::file() << std::endl;
      }

      integratorPtr_->setup();

      // Main MD loop
      timer.start();
      for (int i = 0; i < nStep; ++i) {
         integratorPtr_->step();
      }
      timer.stop();

      // Calculate nAtomTot (correct value only on master).
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
   * Returns total on all processors on master, 0.0 on others.
   */
   double System::kineticEnergy()
   {
      double energy    = 0.0;
      double energyAll = 0.0;
      double mass;
      int typeId;

      // Add kinetic energies of local atoms on this processor
      AtomIterator atomIter;
      atomStorage_.begin(atomIter); 
      for( ; !atomIter.atEnd(); ++atomIter){
         typeId = atomIter->typeId();
         mass   = atomTypes_[typeId].mass();
         energy += mass*(atomIter->velocity().square());
      }
      energy = 0.5*energy;

      #ifdef UTIL_MPI
      // Sum values from all processors.
      domain_.communicator().Reduce(&energy, &energyAll, 1, 
                                    MPI::DOUBLE, MPI::SUM, 0);
      #endif

      return energyAll;
   }

   /*
   * Calculate total potential energy
   * 
   * Returns total on all processors on master, 0.0 on others.
   */
   double System::potentialEnergy() 
   {
      double energy = 0.0;
      energy += pairPotentialEnergy();
      energy += bondPotentialEnergy();
      return energy;
   }

   /*
   * Calculate total nonbonded pair potential energy
   * 
   * Returns total on all processors on master, 0.0 on others.
   */
   double System::pairPotentialEnergy()
   {
      double energy    = 0.0;
      double energyAll = 0.0;

      energy = pairPotential().energy();

      #ifdef UTIL_MPI
      // Sum values from all processors.
      domain_.communicator().Reduce(&energy, &energyAll, 1, 
                                    MPI::DOUBLE, MPI::SUM, 0);
      #endif
      return energyAll;
   }

   /*
   * Calculate total bond potential energy
   * 
   * Returns total on all processors on master, 0.0 on others.
   */
   double System::bondPotentialEnergy()
   {
      double energy = 0.0;
      double energyAll = 0.0;

      energy = bondPotential().energy();

      #ifdef UTIL_MPI
      // Sum values from all processors.
      domain_.communicator().Reduce(&energy, &energyAll, 1, 
                                    MPI::DOUBLE, MPI::SUM, 0);
      #endif
      return energyAll;
   }

   /*
   * Read configuration file on master and distribute atoms.
   *
   * \param filename name of configuration file.
   */
   void System::readConfig(std::string filename)
   {
      assert(configIoPtr_);
      configIoPtr_->readConfig(filename, maskedPairPolicy_);
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

   #if 0
   /*
   * Set pointer to a FileMaster.
   */
   void System::setFileMaster(FileMaster &fileMaster)
   {
      assert(!fileMasterPtr_);
      fileMasterPtr_ = &fileMaster;
   }
   #endif

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
