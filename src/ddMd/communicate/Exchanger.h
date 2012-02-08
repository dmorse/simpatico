#ifndef EXCHANGER_H
#define EXCHANGER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/chemistry/Atom.h>
#include <ddMd/boundary/Boundary.h>
#include <util/containers/FMatrix.h>
#include <util/containers/DPArray.h>

namespace DdMd
{

   class Domain;
   class AtomStorage;
   class BondStorage;
   class Buffer;

   using namespace Util;

   /**
   * Class for exchanging ownership of Atoms between processors.
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
      void associate(const Boundary& boundary, const Domain& domain, 
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
      * Exchange ownership of local atoms just before reneighboring.
      * 
      * This method should be called just before rebuilding the neighbor
      * list on each processor, to exchange ownership of local atoms.
      */
      void exchangeAtoms();

      /**
      * Exchange ghosts.
      * 
      * This method should be called after Exchanger::exchange() 
      * during time steps in which processors exchange atom ownership. 
      * It sends the ghost atom coordinates and atom types, and stores 
      * a record of which atoms were sent in each direction for use in 
      * subsequent calls to communicate().
      */
      void exchangeGhosts();

      /**
      * Update ghost atom coordinates.
      * 
      * This method should be called every time step for which there is
      * no exhange of atom ownership. It communicates ghost coordinates
      * for the same ghosts as those sent by the most recent call to
      * the exchangeGhosts() methods.
      */
      void updateGhosts();

   private:

      /*
      * Arrays of pointers to ghosts to be sent in each direction.
      *
      * The FMatrix(i, j) is a DPArray that stores pointers to the ghost 
      * atoms whose positions are sent during the send along cartesian 
      * axis i (i=0/x, 1/y or 2/z) in direction j (j=0/lower or 1/higher).
      * These 6 arrays are filled by the exchangeGhosts() method and used 
      * by subsequent calls of the updateGhosts() method, until the next
      * call of exchangeGhosts().
      */
      FMatrix< DPArray<Atom>, Dimension, 2>  sendArray_;

      /*
      * Arrays of pointers to ghosts received in each direction.
      *
      * The FMatrix(i, j) is a DPArray that stores pointers to the ghost 
      * atoms that are received during the send along cartesian axis
      * i (i=0/x, 1/y or 2/z) in direction j (j=0/lower or 1/higher).
      * These 6 arrays are filled by the exchangeGhosts() method and used 
      * by subsequent calls of the updateGhosts() method, until the next
      * call of exchangeGhosts().
      */
      FMatrix< DPArray<Atom>, Dimension, 2>  recvArray_;

      // Pointer to associated const Boundary object.
      const Boundary*  boundaryPtr_;

      // Pointer to associated const Domain object.
      const Domain*  domainPtr_;

      // Pointer to associated AtomStorage object.
      AtomStorage* atomStoragePtr_;

      // Pointer to associated AtomStorage object.
      BondStorage* bondStoragePtr_;

      // Pointer to associated buffer object.
      Buffer*  bufferPtr_;

      // Cutoff for pair list (potential cutoff + skin).
      double pairCutoff_;

      // Maximum number of ghosts per transmission.
      int ghostCapacity_;

   };

}
#endif
