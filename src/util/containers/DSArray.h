#ifndef UTIL_DS_ARRAY_H
#define UTIL_DS_ARRAY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/ArrayIterator.h>
#include <util/containers/ConstArrayIterator.h>
#include <util/global.h>

namespace Util
{

   /**
   * Dynamically allocated array with variable logical size.
   *
   * A DSArray < Data > is a wrapper for a dynamically allocated C array,
   * with continuous elements and a logical size that may be less than 
   * or equal to its physical capacity. The logical size is the number 
   * of contiguous elements that have been added using the append() method.
   *
   * \ingroup Array_Module
   */
   template <typename Data>
   class DSArray 
   {

   public:

      /**
      * Constructor.
      */
      DSArray();

      /**
      * Copy constructor.
      *
      *\param other the DSArray to be copied.
      */
      DSArray(const DSArray< Data >& other);
   
      /**
      * Assignment, element by element.
      *
      * Capacity of LHS DSArray must be >= size of RHS DSArray.
      *
      * \param other the RHS DSArray 
      */
      DSArray<Data>& operator=(const DSArray<Data>& other);

      /**
      * Destructor.
      */
      virtual ~DSArray();

      /**
      * Allocates the underlying C array.
      *
      * Throw an exception if the DSArray has already been
      * allocated - A DSArray can only be allocated once.
      *
      * \param capacity number of elements to allocate
      */
      void allocate(int capacity);

      /**
      * Set logical size to zero.
      */
      void clear();

      /**
      * Serialize a DSArray to/from an Archive.
      *
      * \param ar        archive 
      * \param version   archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

      /**
      * Set an ArrayIterator to the beginning of this Array.
      *
      * \param iterator ArrayIterator, initialized on output. 
      */
      void begin(ArrayIterator<Data> &iterator);

      /**
      * Set a ConstArrayIterator to the beginning of this Array.
      *
      * \param iterator ConstArrayIterator, initialized on output. 
      */
      void begin(ConstArrayIterator<Data> &iterator) const;

      /**
      * Mimic C array subscripting.
      *
      * \param  i array index
      * \return reference to element i
      */
      Data& operator[] (int i);

      /**
      * Mimic C array subscripting.
      *
      * \param i array index
      * \return const reference to element i
      */
      const Data& operator[] (int i) const;

      /**
      * Append data to the end of the array.
      *
      * \param data Data to add to end of array.
      */
      void append(const Data &data);

      /**
      * Return physical capacity of array.
      */
      int capacity() const;

      /**
      * Return logical size of this array (i.e., number of elements).
      */
      int size() const;

      /**
      * Return true if the DSArray has been allocated, false otherwise.
      */
      bool isAllocated() const;

   protected:

      /// Array of Data elements.
      Data *data_;

      /// Logical size of array (number of elements used).
      int  size_;

      /// Maxium size of array
      int capacity_;
   };

   // Method definitions

   /*
   * Constructor.
   */
   template <class Data>
   DSArray<Data>::DSArray()
    : data_(0),
      size_(0),
      capacity_(0)
   {}

   /*
   * Copy constructor.
   */
   template <class Data>
   DSArray<Data>::DSArray(const DSArray< Data >& other) 
   {
      // Precondition
      if (!other.isAllocated()) {
         UTIL_THROW("Other DSArray must be allocated.");
       }

      data_     = new Data[other.capacity_];
      size_     = other.size_;
      capacity_ = other.capacity_;
      for (int i = 0; i < size_; ++i) {
         data_[i] = other.data_[i];
      }
   }
   
   /*
   * Assignment, element by element.
   *
   * Capacity of LHS DSArray must be >= size of RHS DSArray.
   */
   template <class Data>
   DSArray<Data>& DSArray<Data>::operator=(const DSArray<Data>& other) 
   {
      // Precondition
      if (!other.isAllocated()) {
         UTIL_THROW("Other DSArray must be allocated.");
      }

      // Check for self assignment
      if (this == &other) return *this;

      if (!isAllocated()) {
         allocate(other.capacity_);
      } else if (capacity_ != other.capacity_) {
         UTIL_THROW("Cannot assign DSArrays of unequal capacity");
      }

      // Copy elements
      for (int i = 0; i < size_; ++i) {
         data_[i] = other[i];
      }
      size = other.size_;

      return *this;
   }

   /*
   * Destructor.
   */
   template <class Data>
   DSArray<Data>::~DSArray()
   {
       if (data_) {
          delete [] data_;
       }
   }

   /*
   * Allocates the underlying C array.
   */
   template <class Data>
   void DSArray<Data>::allocate(int capacity)
   {
      if (data_) {
         UTIL_THROW("Cannot re-allocate a DSArray");
      }
      if (capacity <= 0) {
         UTIL_THROW("Cannot allocate a DSArray with capacity <= 0");
      }
      data_     = new Data[capacity];
      capacity_ = capacity;
   }

   /*
   * Serialize a DSArray to/from an Archive.
   */
   template <class Data>
   template <class Archive>
   void DSArray<Data>::serialize(Archive& ar, const unsigned int version)
   {
      int capacity;
      if (Archive::is_saving()) {
         capacity = capacity_;
      }
      ar & capacity;
      ar & size_;
      if (Archive::is_loading()) {
         if (!isAllocated()) {
            if (capacity > 0) {
               allocate(capacity);
            }
         } else {
            if (size_ > capacity_) {
               UTIL_THROW("Inconsistent DSArray size and capacity");
            }
         }
      }
      for (int i = 0; i < size_; ++i) {
         ar & data_[i];
      }
   }

   /**
   * Set an ArrayIterator to the beginning of this Array.
   *
   * \param iterator ArrayIterator, initialized on output. 
   */
   template <class Data>
   inline
   void DSArray<Data>::begin(ArrayIterator<Data> &iterator)
   {
      iterator.setCurrent(data_);
      iterator.setEnd(data_ + size_);
   }

   /*
   * Set a ConstArrayIterator to the beginning of this Array.
   */
   template <class Data>
   inline 
   void DSArray<Data>::begin(ConstArrayIterator<Data> &iterator) const
   {
      iterator.setCurrent(data_);
      iterator.setEnd(data_ + size_);
   }

   /*
   * Mimic C array subscripting.
   */
   template <class Data>
   inline Data& DSArray<Data>::operator[] (int i)
   {
      assert(i < size_);
      assert(i >= 0);
      return data_[i];
   }

   /*
   * Mimic C array subscripting.
   */
   template <class Data>
   inline const Data& DSArray<Data>::operator[] (int i) const
   {
      assert(i < size_);
      assert(i >= 0 );
      return data_[i];
   }

   /*
   * Append data to the end of the array.
   */
   template <class Data>
   inline void DSArray<Data>::append(const Data &data) 
   {
      if (size_ == capacity_) {
         UTIL_THROW("Attempt to add to full DSArray");
      }
      data_[size_] = data;
      ++size_;
   }

   /*
   * Set logical size to zero.
   */
   template <class Data>
   inline void DSArray<Data>::clear() 
   {  size_ = 0; }

   /*
   * Return physical capacity of array.
   */
   template <class Data>
   inline int DSArray<Data>::capacity() const
   {  return capacity_; }

   /*
   * Return logical size of this array (i.e., number of elements).
   */
   template <class Data>
   inline int DSArray<Data>::size() const
   {  return size_; }

   /*
   * Return true if the DSArray has been allocated, false otherwise.
   */
   template <class Data>
   inline bool DSArray<Data>::isAllocated() const
   {  return !(data_ == 0); }

} 
#endif
