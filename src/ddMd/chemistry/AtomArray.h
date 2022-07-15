#ifndef DDMD_ATOM_ARRAY_H
#define DDMD_ATOM_ARRAY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/Array.h>   // base class template
#include "AtomContext.h"             // context structure.

namespace Util {
   class Vector;
}

namespace DdMd
{

   class Atom;
   class Mask;
   class Plan;

   using namespace Util;

   /**
   * An array of Atom objects.
   *
   * Access to atoms provided by an overloaded subscript ([]) operator,
   * which is inherited from Array<Atom> base class. The interface is
   * similar to that of a DArray, except that copy construction and
   * assignment (=) are prohibited (i.e., private and not implemented).
   *
   * The current implementation places data for some "pseudo-members" of 
   * each Atom in separate arrays. The interface of an Atom hides this, 
   * allowing pseudo-members to be accessed as if they were true class
   * members. See file Atom.h for further implementation details.
   */
   class AtomArray : public Array<Atom>
   { 

   public: 
   
      /**
      * Constructor.
      */
      AtomArray();
   
      /**
      * Destructor.
      */
      virtual ~AtomArray();

      /**
      * Allocate memory on the heap.
      *
      * Throw an Exception if this is already allocated.
      *
      * \param capacity number of elements to allocate.
      */
      void allocate(int capacity); 

      /**
      * Set force vector to zero for all atoms in this array.
      */
      void zeroForces(); 

      /**
      * Return true if this is already allocated, false otherwise.
      */
      bool isAllocated() const;
  
   private:
 
      using Array<Atom>::data_;
      using Array<Atom>::capacity_;

      /*
      * The following C-arrays store data for Atom "psuedo-members".
      * Data associated with the Atom in element i of the main data_
      * array is stored in element i in each of these arrays.
      */

      /**
      * C-array of Atom velocities.
      */
      Vector* velocities_;

      /**
      * C-array of Atom Mask objects.
      */
      Mask* masks_;

      /**
      * C-array of communication Plan data.
      */
      Plan* plans_; 

      /**
      * C-array of unique global atom ids (tags).
      */
      int* ids_;

      /**
      * C-array of bit fields for groups that each Atom belongs to.
      */
      unsigned int* groups_;

      /**
      * C-array of AtomContext structs, one per atom.
      */
      AtomContext* contexts_;

      // Prohibited functions (to suppress compilier generated default)

      /**
      * Copy constructort (prohibited - private and not implemented).
      */
      AtomArray(const AtomArray& other);
   
      /**
      * Assignment (prohibited - private and not implemented).
      */
      AtomArray& operator = (const AtomArray& other); 

   // friends:

      friend class Atom;

   };

}
#endif
