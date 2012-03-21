#ifndef DDMD_ATOM_ARRAY_CPP
#define DDMD_ATOM_ARRAY_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "AtomArray.h"

namespace DdMd
{

   /*
   * Constructor.
   *
   * The data_ and capacity_ are nullified in Array<Atom> constructor.
   */
   AtomArray::AtomArray() 
    : Array<Atom>()
   {}

   /*
   * Destructor.
   */
   AtomArray::~AtomArray()
   {
      if (data_) {
         delete [] data_;
      }
   }

   /*
   * Allocate the underlying C array.
   */
   void AtomArray::allocate(int capacity) 
   {
      if (!(data_ == 0)) {
         UTIL_THROW("Cannot re-allocate an AtomArray");
      }
      if (capacity <= 0) {
         UTIL_THROW("Cannot allocate with capacity <= 0");
      }
      data_     = new Atom[capacity];
      capacity_ = capacity;
   }

   /*
   * Return true if this is already allocated, false otherwise.
   */
   bool AtomArray::isAllocated() const 
   {  return !(data_ == 0); }

}
#endif
