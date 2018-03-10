#ifndef MCMD_MD_SNAPSHOT_H
#define MCMD_MD_SNAPSHOT_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/transition/MdParticle.h>
#include <simp/boundary/Boundary.h>
#include <util/containers/DArray.h>

namespace McMd {

   class MdSystem;

   using namespace Util;
   using namespace Simp;

   class MdSnapShot
   {

      MdSnapShot(MdSystem& system);

      // Copy constructor

      // Assignment

      /**
      * Set state of associated MdSystem to state of this snapshot.
      */
      void setSystemState();

      /**
      * Copy current state of associated MdSystem to this MdSnapShot.
      */
      void getSystemState(int iStep = 0);

      /**
      * Unset state of snapshot, mark as unknown.
      */
      void unset();

      /**
      * Is this snapshot set to a known state?
      */
      bool isSet();

      // Path sampling operations

      /**
      * Reverse all velocities in this snapshot.
      */
      void reverseVelocities();

      /**
      * Add random velocities to a snapshot.
      */
      void addRandomVelocities(double sigma);

      // Accessors

      /**
      * Number of atoms for which memory is allocated;
      */
      int nAtom() const;

      /**
      * Get specific MdParticle by (non-const) reference.
      *
      * \int i particle index
      */
      MdParticle& particle(int i);

      /**
      * Get specific MdParticle by const reference.
      */
      MdParticle const & particle(int i) const;

      /**
      * Get boundary by (non-const) reference.
      */
      Boundary& boundary();

      /**
      * Get boundary by const reference.
      */
      Boundary const & boundary() const;

      /**
      * Get associated system by (non-const) reference.
      */
      MdSystem& system();

      /**
      * Get boundary by const reference.
      */
      MdSystem const & system() const;

      /**
      * Time step of this snapshot.
      */
      int iStep() const;

   private:

      // Array of MdParticles
      DArray<MdParticle> particles_;

      // System Boundary
      Boundary boundary_;

      // Pointer to associated system.
      MdSystem* systemPtr_;

      // Number of atoms in snapshot.
      int nAtom_;

      // Time step stamp
      int iStep_;

      // Is the snapshot state set?
      bool isSet_;

   };

   // Inline Member Functions

   /*
   * Get associated system by reference.
   */
   inline
   MdSystem& MdSnapShot::system()
   {  return *systemPtr_; }

   /*
   * Get associated system by reference.
   */
   inline
   MdSystem const & MdSnapShot::system() const
   {  return *systemPtr_; }

   /*
   * Number of atoms for which memory is allocated;
   */
   inline
   int MdSnapShot::nAtom() const
   {  return nAtom_; }

   /*
   * Get specific MdParticle by (non-const) reference.
   */
   inline
   MdParticle& MdSnapShot::particle(int i)
   {  return particles_[i]; }

   /*
   * Get specific MdParticle by const reference.
   */
   inline
   MdParticle const & MdSnapShot::particle(int i) const
   {  return particles_[i]; }

   /*
   * Get boundary by (non-const) reference.
   */
   inline
   Boundary& MdSnapShot::boundary()
   {  return boundary_; }

   /*
   * Get boundary by const reference.
   */
   inline
   Boundary const & MdSnapShot::boundary() const
   {  return boundary_; }

   /*
   * Time step of this snapshot.
   */
   inline
   int MdSnapShot::iStep() const
   {  return iStep_; }

   /*
   * Unset state of snapshot, mark as unknown.
   */
   inline
   void MdSnapShot::unset()
   {  
      isSet_ = false; 
      nAtom_ = 0; 
   }

   /*
   * Is this snapshot set to a known state?
   */
   inline
   bool MdSnapShot::isSet()
   {  return isSet_; }

}
#endif
