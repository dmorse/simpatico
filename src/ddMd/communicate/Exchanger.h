#ifndef DDMD_EXCHANGER_H
#define DDMD_EXCHANGER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/chemistry/Atom.h>
#include <ddMd/chemistry/Group.h>
#include <util/boundary/Boundary.h>
#include <util/containers/FMatrix.h>
#include <util/containers/APArray.h>

namespace DdMd
{

   class Domain;
   class AtomStorage;
   class BondStorage;
   class Buffer;

   using namespace Util;

   /**
   * Class for exchanging Atoms, Ghosts and Groups between processors.
   *
   * \ingroup DdMd_Communicate_Module
   */
   class Exchanger
   {

   public:

      /**
      * Constructor.
      */
      Exchanger();

      /**
      * Destructor.
      */
      ~Exchanger();

      /**
      * Set pointers to associated objects.
      */
      void associate(const Domain& domain, 
                     const Boundary& boundary, 
                     AtomStorage& atomStorage, 
                     BondStorage& bondStorage, 
                     Buffer& buffer);

      /**
      * Allocate all required memory.
      *
      * Must be called after Buffer is initialized.
      */
      void allocate();

      /**
      * Set width of slab for ghosts.
      *
      * \param pairCutoff cutoff radius for pair list (potential + skin).
      */
      void setPairCutoff(double pairCutoff);

      /**
      * Exchange local atoms and ghosts.
      * 
      * This method should be called just before rebuilding the neighbor
      * list on each processor, to exchange ownership of local atoms and
      * to exchange ghost atoms. The lists of which atoms were sent and
      * received as ghosts by this method are used in subsequent calls
      * to update. 
      */
      void exchange();

      /**
      * Exchange ownership of local atoms.
      *
      * This method exchanges ownershp of local atoms, and calculates
      * communication plan for ghost atoms, but but does not actually
      * exchange ghost atoms.  
      */
      void exchangeAtoms();

      /**
      * Update ghost atom coordinates.
      * 
      * This method should be called every time step for which there is
      * no exhange of atom ownership. It communicates ghost coordinates
      * for the same ghosts as those sent by the most recent call to
      * the exchangeGhosts() methods.
      */
      void update();

   private:

      /**
      * Arrays of pointers to ghosts to be sent in each direction.
      *
      * Element sendArray_(i, j) is a APArray that stores pointers to the
      * ghost atoms whose positions are sent during the send along cartesian 
      * axis i (i=0/x, 1/y or 2/z) in direction j (j=0/lower or 1/higher).
      * These 6 arrays are filled by the exchangeGhosts() method and used 
      * by subsequent calls of the update() method, until the next
      * call of exchangeGhosts().
      */
      FMatrix< APArray<Atom>, Dimension, 2>  sendArray_;

      /**
      * Arrays of pointers to ghosts received in each direction.
      *
      * Element recvArray_(i, j) is a APArray that stores pointers to the 
      * ghost atoms that are received during the send along cartesian axis
      * i (for i=0/x, 1/y or 2/z) in direction j (for j=0/lower or 1/higher).
      * These 6 arrays are filled by the exchangeGhosts() method and used 
      * by subsequent calls of the update() method, until the next
      * call of exchangeGhosts().
      */
      FMatrix< APArray<Atom>, Dimension, 2>  recvArray_;

      /**
      * Array of pointers to empty bonds on this processor.
      * 
      * Used to mark bonds for later removal.
      */
      APArray< Group<2> > emptyBonds_;

      /// Processor boundaries (minima j=0, maxima j=1)
      FMatrix< double, Dimension, 2>  bound_;

      /// Inner boundaries of nonbonded slabs
      FMatrix< double, Dimension, 2>  inner_;

      /// Outer boundaries of nonbonded slabs
      FMatrix< double, Dimension, 2>  outer_;

      /// Pointer to associated const Boundary object.
      const Boundary*  boundaryPtr_;

      /// Pointer to associated const Domain object.
      const Domain*  domainPtr_;

      /// Pointer to associated AtomStorage object.
      AtomStorage* atomStoragePtr_;

      /// Pointer to associated AtomStorage object.
      BondStorage* bondStoragePtr_;

      /// Pointer to associated buffer object.
      Buffer*  bufferPtr_;

      /// Cutoff for pair list (potential cutoff + skin).
      double pairCutoff_;

      /**
      * Exchange ghosts.
      *
      * This method exchanges ghosts, and stores lists of which atoms 
      * are sent and received in each direction for use in subsequent
      * calls to update(). It must be called immediately after
      * exchangeAtoms().
      */
      void exchangeGhosts();

   };

}
#endif
