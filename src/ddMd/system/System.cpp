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
#include <ddMd/interaction/Interaction.h>
#include <ddMd/integrator/Integrator.h>
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
      interactionPtr_(0),
      integratorPtr_(0),
      configIoPtr_(0),
      nAtomType_(0)
   {
      #ifdef UTIL_MPI
      communicatorPtr_ = &communicator;
      domain_.setGridCommunicator(communicator);
      setParamCommunicator(communicator);
      #endif

      interactionPtr_ = new Interaction(*this);
      integratorPtr_  = new Integrator(*this);
      configIoPtr_    = new ConfigIo(*this, buffer_);

      // Set connections between objects
      domain_.setBoundary(boundary_);
      exchanger_.associate(boundary_, domain_, 
                          atomStorage_, bondStorage_, buffer_);

   }

   /**
   * Read parameters, allocate memory and initialize.
   */
   void System::readParam(std::istream& in)
   {

      // Preconditions
      assert(interactionPtr_);
      assert(integratorPtr_);
      assert(configIoPtr_);

      readBegin(in, "System");
      readParamComposite(in, domain_);
      readParamComposite(in, atomStorage_);
      readParamComposite(in, bondStorage_);
      read<int>(in, "nAtomType", nAtomType_);
      pairPotential_.setNAtomType(nAtomType_);
      readParamComposite(in, pairPotential_);
      readParamComposite(in, *interactionPtr_);
      readParamComposite(in, *integratorPtr_);
      readParamComposite(in, random_);
      readParamComposite(in, buffer_);
      readParamComposite(in, *configIoPtr_);

      exchanger_.allocate();
      exchanger_.setPairCutoff(interactionPtr_->cutoff());

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

      integratorPtr_->initialize();

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
         Log::file() << "maxNPair           " << interaction().pairList().maxNPair()
                     << std::endl;
         Log::file() << "maxNAtom           " << interaction().pairList().maxNAtom()
                     << std::endl;
         Log::file() << "buildCounter       " << interaction().pairList().buildCounter()
                     << std::endl;
         Log::file() << "steps / build      "
                     << double(nStep)/double(interaction().pairList().buildCounter())
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

      // Add kinetic energies of local atoms on this processor
      AtomIterator atomIter;
      atomStorage_.begin(atomIter); 
      for( ; !atomIter.atEnd(); ++atomIter){
         energy += 0.5*(atomIter->velocity().square());
      }

      #ifdef UTIL_MPI
      // Sum values from all processors.
      domain_.communicator().Reduce(&energy, &energyAll, 1, 
                                    MPI::DOUBLE, MPI::SUM, 0);
      #endif

      return energyAll;
   }

   /*
   * Calculate total pair potential energy
   * 
   * Returns total on all processors on master, 0.0 on others.
   */
   double System::pairPotentialEnergy()
   {
      double energy    = 0.0;
      double energyAll = 0.0;

      energy = interactionPtr_->pairPotential();

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
         if (sqrt(maxSqDispAll) > 0.5*interactionPtr_->skin()) {
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
      return true; 
   }

}
#endif
