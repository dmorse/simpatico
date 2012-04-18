#ifndef A_P_ARRAY_H
#define A_P_ARRAY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/PArray.h>
#include <util/global.h>

namespace Util
{

   /**
   * Automatically resizable PArray.
   *
   * An APArray is a PArray that grows as needed as objects are appended.
   * Like any PArray, it holds pointers to objects, rather than objects.
   * The associated objects are not destroyed when a PArray is deallocated
   * or destroyed.
   *
   * \ingroup Pointer_Array_Module
   */
   template <typename Data>
   class APArray : public PArray<Data>
   {

   public:

      /**
      * Constructor.
      */
      APArray();

      /**
      * Copy constructor, copy pointers.
      *
      * Allocates new Data* array and copies pointers to Data objects.
      *
      *\param other the APArray to be copied.
      */
      APArray(const APArray<Data>& other);
   
      /**
      * Assignment, element by element.
      *
      * Preconditions: 
      * - Both this and other APArrays must be allocated.
      * - Capacity of this APArray must be >= size of RHS APArray.
      *
      * \param other the rhs APArray 
      */
      APArray<Data>& operator=(const APArray<Data>& other);

      /**
      * Destructor.
      *
      * Deletes array of pointers, if allocated previously.
      * Does not delete the associated Data objects.
      */
      virtual ~APArray();

      /**
      * Append an element to the end of the sequence.
      *
      * Resizes array if space is inadequate. 
      *
      * \param data Data object to be appended
      */
      void append(Data& data);

      /**
      * Reserve memory for specified number of elements.
      *
      * Resizes and copies array if requested capacity is less than the
      * current capacity. Does nothing if requested capacity is greater
      * than current capacity.
      *
      * \param capacity number of elements for which to reserve space.
      */
      void reserve(int capacity);

      /**
      * Deallocate (delete) underlying array of pointers.
      */
      void deallocate();

      /**
      * Reset to empty state.
      */ 
      void clear();

   protected:

      using PArray<Data>::ptrs_;
      using PArray<Data>::capacity_;
      using PArray<Data>::size_;

   };

   /*
   * Default constructor.
   */
   template <typename Data>
   inline APArray<Data>::APArray()
    : PArray<Data>()
   {}

   /**
   * Copy constructor, copy pointers.
   *
   * Allocates a new Data* array and copies all pointer values.
   *
   *\param other the APArray to be copied.
   */
   template <typename Data>
   APArray<Data>::APArray(const APArray<Data>& other) 
    : PArray<Data>()
   {
      if (other.ptrs_ == 0) {

         ptrs_ = 0;
         capacity_ = 0;
         size_ = 0;

      } else { 

         // Allocate array of Data* pointers
         ptrs_  = new Data*[other.capacity_];
         capacity_ = other.capacity_;
         size_ = other.size_;

         // Copy pointers
         int i;
         for (i = 0; i < size_; ++i) {
            ptrs_[i] = other.ptrs_[i];
         }

         // Nullify unused elements of ptrs_ array
         if (capacity_ > size_) {
            for (i = size_; i < capacity_; ++i) {
               ptrs_[i] = 0;
            }
         }

      }
   }

   /*
   * Assignment, element by element.
   */
   template <typename Data>
   APArray<Data>& APArray<Data>::operator=(const APArray<Data>& other) 
   {
      // Check for self assignment
      if (this == &other) return *this;

      clear();
      for (int i = 0; i < other.size_; ++i) {
         append(other[i]);
      }
      return *this;
   }

   // Destructor.
   template <typename Data>
   inline APArray<Data>::~APArray()
   {
      if (ptrs_) {
         delete [] ptrs_;
      }
   }

   /*
   * Reserve space for the underlying array of Data* pointers.
   */
   template <typename Data>
   void APArray<Data>::reserve(int capacity) 
   {
      if (capacity <= 0) {
         UTIL_THROW("Cannot reserve with capacity <=0");
      }
      if (ptrs_ == 0) {
         ptrs_ = new Data*[capacity];
         capacity_ = capacity;
      } else if (capacity > capacity_) {
         Data** old = ptrs_;
         ptrs_ = new Data*[capacity];
         capacity_ = capacity;
         for (int i = 0; i < size_; ++i) {
            ptrs_[i] = old[i];
         }
         delete [] old;
      }
   }

   /*
   * Deallocate associated memory.
   */
   template <typename Data>
   void APArray<Data>::deallocate() 
   {  
      if (ptrs_) {
         delete [] ptrs_;
         capacity_ = 0; 
         size_ = 0; 
      }
   }

   /*
   * Append an element to the end of the PArray.
   */
   template <typename Data>
   void APArray<Data>::append(Data& data) 
   {
      if (ptrs_ == 0) {
         ptrs_ = new Data*[64];
         capacity_ = 64;
         size_ = 0;
      } else if (size_ == capacity_) {
         Data** old = ptrs_;
         capacity_ = 2*capacity_;
         ptrs_ = new Data*[capacity_];
         for (int i = 0; i < size_; ++i) {
            ptrs_[i] = old[i];
         }
         delete [] old;
      }
      ptrs_[size_] = &data;
      ++size_;
   }

   /*
   * Reset to empty state, without deallocating.
   */
   template <typename Data>
   void APArray<Data>::clear()
   {
      size_ = 0;
   }



} 
#endif
