#ifndef DDMD_ATOM_ARRAY_H
#define DDMD_ATOM_ARRAY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

//#include "Atom.h"
#include <util/containers/Array.h>

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
   * The current implementation of this class is almost identical to
   * DArray<Atom>, except for a prohibition on copying and assignment.
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
      * Return true if this is already allocated, false otherwise.
      */
      bool isAllocated() const;
  
   private:
 
      using Array<Atom>::data_;
      using Array<Atom>::capacity_;

      Vector *velocities_;
      Mask   *masks_;
      Plan   *plans_; 
      int    *ids_;

      /**
      * Copy ctor (prohibited - private and not implemented).
      */
      AtomArray(const AtomArray& other);
   
      /**
      * Assignment (prohibited - private and not implemented).
      */
      AtomArray& operator = (const AtomArray& other); 

      friend class Atom;

   }; 

}
#endif
