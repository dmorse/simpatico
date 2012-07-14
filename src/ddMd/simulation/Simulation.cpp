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
#include <ddMd/integrators/Integrator.h>
#include <ddMd/integrators/IntegratorFactory.h>
#include <ddMd/configIos/ConfigIo.h>
#include <ddMd/configIos/ConfigIoFactory.h>
#include <ddMd/configIos/DdMdConfigIo.h>
#include <ddMd/util/FileMaster.h>
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
#include <ddMd/potentials/angle/AnglePotentialImpl.h>
#include <ddMd/potentials/angle/AngleFactory.h>
#endif
#ifdef INTER_DIHEDRAL
#include <ddMd/potentials/dihedral/DihedralPotential.h>
#include <ddMd/potentials/dihedral/DihedralPotentialImpl.h>
#include <ddMd/potentials/dihedral/DihedralFactory.h>
#endif
#ifdef INTER_EXTERNAL
#include <ddMd/potentials/external/ExternalPotential.h>
#include <ddMd/potentials/external/ExternalPotentialImpl.h>
#include <ddMd/potentials/external/ExternalFactory.h>
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
#include <unistd.h>

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
      reverseUpdateFlag_(false)
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
                           atomStorage_, bondStorage_, 
                           #ifdef INTER_ANGLE
                           angleStorage_,
                           #endif
                           #ifdef INTER_DIHEDRAL
                           dihedralStorage_,
                           #endif
                           buffer_);

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
      //if (logFile_.is_open()) logFile_.close();
      #endif
   }

   /*
   * Process command line options.
   */
   void Simulation::setOptions(int argc, char **argv)
   {
      bool  eflag  = false;
   
      // Read command-line arguments
      int c;
      opterr = 0;
      while ((c = getopt(argc, argv, "epr:")) != -1) {
         switch (c) {
         case 'e':
           eflag = true;
           break;
         case '?':
           std::cout << "Unknown option -" << optopt << std::endl;
         }
      }
   
      if (eflag) {
         // Enable echoing of parameters to log file as they are read.
         Util::ParamComponent::setEcho(true);
      }
   }

   /**
   * Read parameters, allocate memory and initialize.
   */
   void Simulation::readParam(std::istream& in)
   {
      // Preconditions
      assert(pairPotentialPtr_ == 0);
      assert(bondPotentialPtr_ == 0);
      assert(integratorPtr_ == 0);
      assert(configIoPtr_ == 0);

      readBegin(in, "Simulation");

      readParamComposite(in, domain_);

      readFileMaster(in);

      // Read types
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

      #ifndef DDMD_NOPAIR
      // Pair Potential
      pairPotentialPtr_ = pairFactory().factory(pairStyle());
      pairPotentialPtr_->setNAtomType(nAtomType_);
      readParamComposite(in, *pairPotentialPtr_);
      pairPotentialPtr_->setForceCommFlag(reverseUpdateFlag_);
      #endif

      // Bond Potential
      bondPotentialPtr_ = bondFactory().factory(bondStyle());
      bondPotentialPtr_->setNBondType(nBondType_);
      readParamComposite(in, *bondPotentialPtr_);

      #ifdef INTER_ANGLE
      // Angle potential
      if (nAngleType_) {
         anglePotentialPtr_ = angleFactory().factory(angleStyle());
         anglePotentialPtr_->setNAngleType(nAngleType_);
         readParamComposite(in, *anglePotentialPtr_);
      }
      #endif

      #ifdef INTER_DIHEDRAL
      // Dihedral potential
      if (nDihedralType_) {
         dihedralPotentialPtr_ = dihedralFactory().factory(dihedralStyle());
         dihedralPotentialPtr_->setNDihedralType(nDihedralType_);
         readParamComposite(in, *dihedralPotentialPtr_);
      }
      #endif

      #ifdef INTER_EXTERNAL
      // External potential
      if (hasExternal_) {
         externalPotentialPtr_ = externalFactory().factory(externalStyle());
         externalPotentialPtr_->setNAtomType(nAtomType_);
         readParamComposite(in, *externalPotentialPtr_);
      }
      #endif

      readEnsembles(in);

      // Integrator
      std::string className;
      bool isEnd;
      integratorPtr_ = 
         integratorFactory().readObject(in, *this, className, isEnd);
      if (!integratorPtr_) {
         std::string msg("Unknown Integrator subclass name: ");
         msg += className;
         UTIL_THROW("msg.c_str()");
      }

      readParamComposite(in, random_);
      readParamComposite(in, *diagnosticManagerPtr_);

      exchanger_.setPairCutoff(pairPotentialPtr_->cutoff());
      exchanger_.allocate();

      readEnd(in);
   }

   /**
   * If no FileMaster exists, create and initialize one. 
   */
   void Simulation::readFileMaster(std::istream &in)
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
   * Read EnergyEnsemble and BoundaryEnsemble
   */
   void Simulation::readEnsembles(std::istream &in)
   {
      readParamComposite(in, *energyEnsemblePtr_);
      readParamComposite(in, *boundaryEnsemblePtr_);
   }

   /*
   * Read and implement commands in an input script.
   */
   void Simulation::readCommands(std::istream &in)
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
         } else
         if (command == "THERMALIZE") {
            double temperature;
            inBuffer >> temperature;
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
         } else
         if (command == "WRITE_PARAM") {
            inBuffer >> filename;
            fileMaster().openOutputFile(filename, outputFile);
            writeParam(outputFile);
            outputFile.close();
         } else
         if (command == "SET_CONFIG_IO") {
            std::string classname;
            inBuffer >> classname;
            setConfigIo(classname);
         } else
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
   void Simulation::readCommands()
   {  readCommands(fileMaster().commandFile()); }

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
      pairPotential().addForces();
      bondPotential().addForces();
      #ifdef INTER_ANGLE
      if (nAngleType_) {
         anglePotential().addForces();
      }
      #endif
      #ifdef INTER_DIHEDRAL
      if (nDihedralType_) {
         dihedralPotential().addForces();
      }
      #endif
      #ifdef INTER_EXTERNAL
      if (hasExternal_) {
         externalPotential().addForces();
      }
      #endif

      // Reverse communication (if any)
      if (reverseUpdateFlag_) {
         exchanger_.reverseUpdate();
      }

   }

   /*
   * Integrate.
   */
   void Simulation::simulate(int nStep)
   {
      // Preconditions
      assert(integratorPtr_);

      if (domain_.isMaster()) {
         Log::file() << std::endl;
      }

      integratorPtr_->setup();
      integratorPtr_->run(nStep);
      integratorPtr_->outputStatistics(Log::file());
   }

   /*
   * Calculate total kinetic energy
   *
   * Call on all processors. 
   */
   void Simulation::computeKineticEnergy()
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
   double Simulation::kineticEnergy()
   {  return kineticEnergy_; }

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
   * Compute all potential energy contributions.
   */
   double Simulation::potentialEnergy() 
   {
      double energy = 0.0;
      energy += pairPotential().energy();
      energy += bondPotential().energy();
      #ifdef INTER_ANGLE
      if (nAngleType_) {
         energy += anglePotential().energy();
      }
      #endif
      #ifdef INTER_DIHEDRAL
      if (nDihedralType_) {
         energy += dihedralPotential().energy();
      }
      #endif
      #ifdef INTER_EXTERNAL
      if (hasExternal_) {
         energy += externalPotential().energy();
      }
      #endif
      return energy;
   }

   /*
   * Read configuration file on master and distribute atoms.
   */
   void Simulation::readConfig(const std::string& filename)
   {
      std::ifstream inputFile;
      if (domain_.isMaster()) {
         fileMaster().openInputFile(filename, inputFile);
      }
      if (configIoPtr_ == 0) {
         configIoPtr_ = new DdMdConfigIo(*this);
         configIoPtr_->associate(domain_, boundary_,
                                 atomStorage_, bondStorage_, 
                                 #ifdef INTER_ANGLE
                                 angleStorage_,
                                 #endif
                                 #ifdef INTER_DIHEDRAL
                                 dihedralStorage_,
                                 #endif
                                 buffer_);
         configIoPtr_->initialize();
      }
      configIoPtr_->readConfig(inputFile, maskedPairPolicy_);
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
      if (configIoPtr_ == 0) {
         configIoPtr_ = new DdMdConfigIo(*this);
         configIoPtr_->associate(domain_, boundary_,
                                 atomStorage_, bondStorage_, 
                                 #ifdef INTER_ANGLE
                                 angleStorage_,
                                 #endif
                                 #ifdef INTER_DIHEDRAL
                                 dihedralStorage_,
                                 #endif
                                 buffer_);
         configIoPtr_->initialize();
      }
      configIoPtr_->writeConfig(outputFile);
      if (domain_.isMaster()) {
         outputFile.close();
      }
   }

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

   // ConfigIoIo Management

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
      configIoPtr_->associate(domain_, boundary_,
                              atomStorage_, bondStorage_, 
                              #ifdef INTER_ANGLE
                              angleStorage_,
                              #endif
                              #ifdef INTER_DIHEDRAL
                              dihedralStorage_,
                              #endif
                              buffer_);
      configIoPtr_->initialize();
   }

   /**
   * Return true if this Simulation is valid, or throw an Exception.
   */
   bool Simulation::isValid()
   {
      atomStorage_.isValid();

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

      #else

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

      return true; 
   }

   /*
   * Set flag to specify if reverse communication is enabled.
   */
   void Simulation::setForceCommFlag(bool reverseUpdateFlag)
   {  
      reverseUpdateFlag_ = reverseUpdateFlag; 
      if (pairPotentialPtr_) {
         pairPotentialPtr_->setForceCommFlag(reverseUpdateFlag);
      }
   }

}
#endif
