/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/transition/MdSnapShot.h>
#include <mcMd/mdSimulation/MdSystem.h>
#include <mcMd/mdSimulation/MdSimulation.h>
#include <util/global.h>

namespace McMd {

   using namespace Util;

   MdSnapShot::MdSnapShot(MdSystem& system)
    : particles_(),
      boundary_(),
      systemPtr_(&system),
      iStep_(0),
      isSet_(false)
   {}

   // Copy constructor

   // Assignment

   /*
   * Store current state of system in this MdSnapShot.
   */
   void MdSnapShot::getSystemState(int iStep)
   {
      int nAtom_ = system().nAtom();
      UTIL_CHECK(nAtom_ > 0);
      if (particles_.capacity() == 0) {
         particles_.allocate(nAtom_);
      } else {
         UTIL_CHECK(nAtom_ == particles_.capacity());
      }

      // Simulation& simulation = system().simulation();
      System::MoleculeIterator molIter;
      Molecule::AtomIterator atomIter;
      int iSpecies, iAtom;
      int nSpecies = system().simulation().nSpecies();

      // 1st half velocity Verlet, loop over atoms 
      iAtom = 0;
      for (iSpecies=0; iSpecies < nSpecies; ++iSpecies) {
         system().begin(iSpecies, molIter); 
         for ( ; molIter.notEnd(); ++molIter) {
            for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
               particles_[iAtom].position = atomIter->position();
               particles_[iAtom].velocity = atomIter->velocity();
               ++iAtom;
            }
         }
      }
      UTIL_CHECK(iAtom == nAtom_);

      iStep_ = iStep;
      isSet_ = true;
   }

   /*
   * Set system to state of this snapshot.
   */
   void MdSnapShot::setSystemState()
   {
      UTIL_CHECK(nAtom_ > 0);
      UTIL_CHECK(nAtom_ == system().nAtom());
      UTIL_CHECK(nAtom_ == particles_.capacity());

      // Simulation& simulation = system().simulation();
      System::MoleculeIterator molIter;
      Molecule::AtomIterator atomIter;
      int iSpecies, iAtom;
      int nSpecies = system().simulation().nSpecies();

      // 1st half velocity Verlet, loop over atoms 
      iAtom = 0;
      for (iSpecies=0; iSpecies < nSpecies; ++iSpecies) {
         system().begin(iSpecies, molIter); 
         for ( ; molIter.notEnd(); ++molIter) {
            for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
               atomIter->position() = particles_[iAtom].position;
               atomIter->velocity() = particles_[iAtom].velocity;
               ++iAtom;
            }
         }
      }
      UTIL_CHECK(iAtom == nAtom_);
   }

   // Path sampling operations
   /*
   * Reverse all velocities in this snapshot.
   */
   void MdSnapShot::reverseVelocities()
   {
      UTIL_CHECK(isSet_);
      UTIL_CHECK(nAtom_ > 0);
      for (int i = 0; i < nAtom_; ++i) {
         particles_[i].velocity *= -1;
      }
   }

   /*
   * Add random velocities to a snapshot.
   */
   void MdSnapShot::addRandomVelocities(double sigma)
   {
   }

}
