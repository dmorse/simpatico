#ifndef UTIL_MEMORY_H
#define UTIL_MEMORY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>
#include <stddef.h>
#include <iostream>
#include <new>

namespace Util
{

   /**
   * Provides method to allocate array.
   *
   * The Memory::allocate() method invokes the new operator within
   * a try catch block, and keeps track of the total memory 
   * allocated.
   */
   class Memory
   { 
   public:

      /**
      * Allocate a C array.
      *
      * Allocates a Data array of size elements, assigns ptr the
      * address of the first element. 
      */
      template <typename Data>
      static void allocate(Data*& ptr, size_t size);

      /**
      * Allocate a C array.
      *
      * Allocates a Data array of size elements, assigns ptr the
      * address of the first element. 
      */
      template <typename Data>
      static void deallocate(Data*& ptr, size_t size);

      /**
      * Return total amount of memory currently allocated.
      */
      static size_t total();

      /**
      * Return the maximum amount of allocated heap memory thus far.
      *
      * This function returns the temporal maximum of total().
      */
      static size_t max();

      #ifdef UTIL_MPI
      /**
      * Return max for any processor in communicator.
      */
      static size_t max(MPI::Intracomm& communicator);
      #endif

      /**
      * Call this just to guarantee initialization of static memory.
      */
      static void initStatic();
   
   private: 

      /// Total amount of memory allocated, in bytes. 
      static size_t total_;
   
      /// Maximum amount of memory allocated, in bytes. 
      static size_t max_;
   
   };
   
   /*
   * Allocate a C array.
   */
   template <typename Data>
   void Memory::allocate(Data*& ptr, size_t size)
   {
      try {
         ptr  = new Data[size];
      } catch (std::bad_alloc&) {
         std::cout << "Allocation error" << std::endl;
         throw;
      }
      total_ += size*sizeof(Data);
      if (total_ > max_) max_ = total_;
   }

   /*
   * De-allocate a C array.
   */
   template <typename Data>
   void Memory::deallocate(Data*& ptr, size_t size)
   {
      if (ptr) {
         delete [] ptr;
         total_ -= size*sizeof(Data);
         ptr = 0;
      }
   }

} 
#endif
