#ifndef UTIL_D_ARRAY_H
#define UTIL_D_ARRAY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/Array.h>
#include <util/global.h>


namespace Util
{

   /**
   * Dynamically allocatable Array class template.
   *
   * A DArray wraps a dynamically allocated C Array. Unlike an STL 
   * std::vector, a DArray cannot be resized after it is allocated.
   *
   * \ingroup Array_Module
   */
   template <typename Data>
   class DArray : public Array<Data>
   { 

      using Array<Data>::data_;
      using Array<Data>::capacity_;
   
   public: 
   
      /**
      * Default constructor.
      */
      DArray();
   
      /**
      * Copy constructor.
      *
      * Allocates new memory and copies all elements by value.
      *
      *\param other the DArray to be copied.
      */
      DArray(const DArray<Data>& other);
   
      /**
      * Assignment operator.
      *
      * If this DArray is not allocated, allocates and copies all elements.
      *
      * If this and the other DArray are both allocated, the capacities must
      * be exactly equal. If so, this method copies all elements.
      *
      * \param other the RHS DArray
      */
      DArray<Data>& operator = (const DArray<Data>& other); 

      /**
      * Destructor.
      *
      * Deletes underlying C array, if allocated previously.
      */
      virtual ~DArray();

      /**
      * Allocate the underlying C array.
      *
      * Throw an Exception if the DArray is already allocated.
      *
      * \param capacity number of elements to allocate.
      */
      void allocate(int capacity); 

      /**
      * Dellocate the underlying C array.
      *
      * Throw an Exception if the DArray is not allocated.
      */
      void deallocate(); 

      /**
      * Return true if the DArray has been allocated, false otherwise.
      */
      bool isAllocated() const;

      /*
      * Serialize a DArray to/from an Archive.
      *
      * \param ar        archive 
      * \param container dynamic array 
      * \param version   archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

   }; 


   /*
   * Constructor.
   */
   template <class Data>
   DArray<Data>::DArray() 
    : Array<Data>()
   {}

   /*
   * Copy constructor.
   *
   * Allocates new memory and copies all elements by value.
   *
   *\param other the DArray to be copied.
   */
   template <class Data>
   DArray<Data>::DArray(const DArray<Data>& other) 
    : Array<Data>()
   {
      // Precondition
      if (!other.isAllocated()) {
         UTIL_THROW("Other DArray must be allocated.");
      }

      data_     = new Data[other.capacity_];
      capacity_ = other.capacity_;
      for (int i = 0; i < capacity_; ++i) {
         data_[i] = other.data_[i];
      }

   }

   /*
   * Destructor.
   */
   template <class Data>
   DArray<Data>::~DArray()
   {
      if (data_) {
         delete [] data_;
      }
   }

   /*
   * Assignment, element-by-element.
   *
   * This operator will allocate memory if not allocated previously.
   *
   * \throw Exception if other DArray is not allocated.
   * \throw Exception if both DArrays are allocated with unequal capacities.
   *
   * \param other the rhs DArray 
   */
   template <class Data>
   DArray<Data>& DArray<Data>::operator = (const DArray<Data>& other) 
   {
      // Precondition
      if (!other.isAllocated()) {
         UTIL_THROW("Other DArray must be allocated.");
      }

      // Check for self assignment
      if (this == &other) return *this;

      if (!isAllocated()) {
         allocate(other.capacity());
      } else if (capacity_ != other.capacity_) {
         UTIL_THROW("Cannot assign DArrays of unequal capacity");
      }

      // Copy elements
      for (int i = 0; i < capacity_; ++i) {
         data_[i] = other[i];
      }

      return *this;
   }

   /*
   * Allocate the underlying C array.
   *
   * Throw an Exception if the DArray has already allocated.
   *
   * \param capacity number of elements to allocate.
   */
   template <class Data>
   void DArray<Data>::allocate(int capacity) 
   {
      if (!(data_ == 0)) {
         UTIL_THROW("Cannot re-allocate a DArray");
      }
      if (capacity <= 0) {
         UTIL_THROW("Cannot allocate a DArray with capacity <= 0");
      }
      data_     = new Data[capacity];
      capacity_ = capacity;
   }

   /*
   * Deallocate the underlying C array.
   *
   * Throw an Exception if this DArray is not allocated.
   */
   template <class Data>
   void DArray<Data>::deallocate() 
   {
      if (!data_) {
         UTIL_THROW("Array is not allocated");
      }
      delete [] data_;
      data_ = 0;
      capacity_ = 0;
   }

   /*
   * Return true if the DArray has been allocated, false otherwise.
   */
   template <class Data>
   inline bool DArray<Data>::isAllocated() const 
   {  return !(data_ == 0); }

   /*
   * Serialize a DArray to/from an Archive.
   */
   template <class Data>
   template <class Archive>
   void DArray<Data>::serialize(Archive& ar, const unsigned int version)
   {
      int capacity;
      if (Archive::is_saving()) {
         capacity = capacity_;
      }
      ar & capacity;
      if (Archive::is_loading()) {
         if (!isAllocated()) {
            if (capacity > 0) {
               allocate(capacity);
            }
         } else {
            if (capacity != capacity_) {
               UTIL_THROW("Inconsistent DArray capacities");
            }
         }
      }
      for (int i = 0; i < capacity_; ++i) {
         ar & data_[i];
      }
   }

} 
#endif
