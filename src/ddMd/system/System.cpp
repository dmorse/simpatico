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
#include <ddMd/potentials/pair/PairInteraction.h>
#include <ddMd/potentials/bond/BondPotential.h>
#include <ddMd/potentials/bond/BondPotentialImpl.h>
#include <ddMd/potentials/bond/BondInteraction.h>
#include <ddMd/integrator/NveIntegrator.h>
#include <ddMd/configIo/ConfigIo.h>
#include <util/util/Log.h>
#include <util/util/Timer.h>

//#include <util/format/Int.h>
//#include <util/format/Dbl.h>
//#include <util/format/Str.h>
//#include <util/util/ioUtil.h>


namespace DdMd
{

   using namespace Util;

   /*
   * Default constructor.
   */
   #ifdef UTIL_MPI
   System::System(MPI::Intracomm& communicator)
    : communicatorPtr_(0),
   #else
   System::System() :
   #endif
      pairPotentialPtr_(0),
      bondPotentialPtr_(0),
      integratorPtr_(0),
      configIoPtr_(0),
      nAtomType_(0),
      nBondType_(0)
   {
      #ifdef UTIL_MPI
      communicatorPtr_ = &communicator;
      domain_.setGridCommunicator(communicator);
      setParamCommunicator(communicator);
      #endif

      // Set connections between objects
      domain_.setBoundary(boundary_);
      exchanger_.associate(domain_, boundary_,
                           atomStorage_, bondStorage_, buffer_);

      // Note: The following objects will become polymorphic,
      // and will be created in readParam by factories.

      pairPotentialPtr_ = new PairPotentialImpl<PairInteraction>(*this);
      bondPotentialPtr_ = new BondPotentialImpl<BondInteraction>(*this);
      integratorPtr_  = new NveIntegrator(*this);
      configIoPtr_  = new ConfigIo();
      configIoPtr_->associate(domain_, boundary_,
                              atomStorage_, bondStorage_, buffer_);

   }

   /**
   * Read parameters, allocate memory and initialize.
   */
   void System::readParam(std::istream& in)
   {

      // Preconditions
      assert(pairPotentialPtr_);
      assert(bondPotentialPtr_);
      assert(integratorPtr_);
      assert(configIoPtr_);

      readBegin(in, "System");
      readParamComposite(in, domain_);
      read<int>(in, "nAtomType", nAtomType_);
      read<int>(in, "nBondType", nBondType_);
      atomTypes_.allocate(nAtomType_);
      for (int i = 0; i < nAtomType_; ++i) {
         atomTypes_[i].setId(i);
      }
      readDArray<AtomType>(in, "atomTypes", atomTypes_, nAtomType_);
      readParamComposite(in, atomStorage_);
      readParamComposite(in, bondStorage_);
      pairPotentialPtr_->setNAtomType(nAtomType_);
      bondPotentialPtr_->setNBondType(nBondType_);
      //readParamComposite(in, pairInteraction_);
      //bondInteraction_.setNBondType(nBondType_);
      //readParamComposite(in, bondInteraction_);
      readParamComposite(in, *pairPotentialPtr_);
      readParamComposite(in, *bondPotentialPtr_);
      readParamComposite(in, *integratorPtr_);
      readParamComposite(in, random_);
      readParamComposite(in, buffer_);
      readParamComposite(in, *configIoPtr_);

      exchanger_.allocate();
      exchanger_.setPairCutoff(pairPotentialPtr_->cutoff());

      readEnd(in);
   }

   /*
   * Set forces on all local atoms to zero.
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
      int nAtomTot = nAtomTotal();

      if (isMaster) {

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
      configIoPtr_->readConfig(filename);
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

   /**
   * Return total number of atoms on all processors.
   */
   int System::nAtomTotal() const
   {
      int nAtom = atomStorage_.nAtom();
      int nAtomAll;
      #ifdef UTIL_MPI
      domain_.communicator().Reduce(&nAtom, &nAtomAll, 1, 
                                    MPI::INT, MPI::SUM, 0);
      #else
      nAtomAll = nAtom;
      #endif
      return nAtomAll;
   }

   /**
   * Return total number of ghosts on all processors.
   * 
   * Reduce operation: Must be called on all processors but returns
   * correct value only on processor 0 of grid communicator.
   */
   int System::nGhostTotal() const
   {
      int nGhost = atomStorage_.nGhost();
      int nGhostAll = 0;
      #ifdef UTIL_MPI
      domain_.communicator().Reduce(&nGhost, &nGhostAll, 1, 
                                    MPI::INT, MPI::SUM, 0);
      #else
      nGhostAll = nGhost;
      #endif
      return nGhostAll;
   }

   /**
   * Return true if this System is valid, or throw an Exception.
   */
   bool System::isValid()
   {
      atomStorage_.isValid();

      // Determine if there are any ghosts on any processor
      int nGhost = nGhostTotal();
      #if UTIL_MPI
      bcast(domain_.communicator(), nGhost, 0);
      #endif
      bool hasGhosts = bool(nGhost);

      bondStorage_.isValid(atomStorage_, domain_.communicator(), hasGhosts);
      return true; 

   }

}
#endif
