#ifndef MCMD_MD_SNAPSHOT_H
#ifndef MCMD_MD_SNAPSHOT_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

namespace McMd {

   class MdSnapShot 
   {

      MdSnapShot(MdSystem& );

      // Copy constructor

      // Assignment

      allocate(nAtom);

      /**
      * Set current state of system to state of this snapshot.
      */
      setSystemState();

      /**
      * Set snapshot state to current state of system.
      */
      getSystemState(int iStep = 0);

      /**
      * Unset state of snapshot, mark as unknown.
      */
      unset()

      /**
      * Is this snapshot set to a state?
      */
      bool isSet();

      // Path sampling operations

      /**
      * Add random velocities to a snapshot.
      */ 
      void addRandomVelocities(double sigma, MdSnapShot const & in);

      /**
      * Copy a snapshot with reversed velocities. (Needed)?
      */ 
      void reverseVelocities(MdSnapShot const & in);

      // Accessors

      /**
      * Number of atoms for which memory is allocated;
      */
      int nAtom() const;

      /**
      * Get specific position const reference.
      *
      * \int 
      */
      MdParticle& particle(int i);

      /**
      * Get specific atom by const reference.
      */
      MdParticle const & particle(int i);
      
      /**
      * Get boundary by reference.
      */ 
      Boundary& boundary();

      /**
      * Get boundary by const reference.
      */ 
      Boundary const & boundary() const;

      /**
      * Time step of this snapshot.
      */
      int iStep() const;

   private:

      // Array of MdParticles
      DArray<MdParticle> particles_;

      // System Boundary
      Boundary boundary_;

      // Time step stamp
      int iStep_;

   };

}
#endif
