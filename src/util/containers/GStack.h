#ifndef UTIL_G_STACK_H
#define UTIL_G_STACK_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/GStack.h>
#include <util/global.h>

namespace Util
{

   /**
   * An automatically growable Stack.
   *
   * A GStack is stack that is implemented as a growable pointer array.
   * Like any pointer array it holds pointers to objects, rather than 
   * objects, and associated objects are not destroyed when the container
   * is deallocated or destroyed.
   *
   * \ingroup Pointer_Array_Module
   */
   template <typename Data>
   class GStack 
   {

   public:

      /**
      * Constructor.
      */
      GStack();

      /**
      * Copy constructor, copy pointers.
      *
      * Allocates new Data* array and copies pointers to Data objects.
      *
      *\param other the GStack to be copied.
      */
      GStack(const GStack<Data>& other);
   
      /**
      * Destructor.
      *
      * Deletes array of pointers, if allocated previously.
      * Does not delete the associated Data objects.
      */
      ~GStack();

      /**
      * Assignment, element by element.
      *
      * Preconditions: 
      * - Both this and other GStacks must be allocated.
      * - Capacity of this GStack must be >= size of RHS GStack.
      *
      * \param other the rhs GStack 
      */
      GStack<Data>& operator=(const GStack<Data>& other);

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

      /**
      * Push an element onto the stack.
      *
      * Resizes array if space is inadequate. 
      *
      * \param data element to be added to stack.
      */
      void push(Data& data);

      /**
      * Pop an element off the stack.
      * 
      * Returns the top element by reference and removes it,
      * decrementing the size by one.
      *
      * \return the top element (which is popped off stack).
      */
      Data& pop();

      /**
      * Return a reference to the top element (don't pop).
      */
      Data& peek();

      /**
      * Return a const ref to the top element (don't pop).
      */
      const Data& peek() const;

      /**
      * Return allocated size.
      *
      * \return Number of elements allocated in array.
      */
      int capacity() const;

      /**
      * Return logical size.
      *
      * \return logical size of this array.
      */
      int size() const;

      /**
      * Is this GStack allocated?
      */
      bool isAllocated() const;

      /**
      * Is this GStack in a valid internal state?
      */
      bool isValid() const;

   private:

      /// GStack of of pointers to Data objects.
      Data** ptrs_;

      /// Allocated size of ptrs_ array.
      int capacity_;

      /// Logical size (number of elements with initialized data).
      int size_;

   };

   /*
   * Constructor.
   */
   template <typename Data>
   inline GStack<Data>::GStack()
    : ptrs_(0),
      capacity_(0),
      size_(0)
   {}

   /**
   * Copy constructor, copy pointers.
   *
   * Allocates a new Data* array and copies all pointer values.
   *
   *\param other the GStack to be copied.
   */
   template <typename Data>
   GStack<Data>::GStack(const GStack<Data>& other) 
    : ptrs_(0),
      capacity_(0),
      size_(0)
   {
      assert(other.capacity_ >= other.size_);
      if (other.ptrs_ == 0) {
         assert(other.capacity_ == 0);
         assert(other.size_ == 0);
         ptrs_ = 0;
         capacity_ = 0;
         size_ = 0;
      } else { 
         assert(other.capacity_ > 0);
         // Allocate array of Data* pointers
         Memory::allocate<Data*>(ptrs_, other.capacity_);
         capacity_ = other.capacity_;
         size_ = other.size_;
         // Copy pointers
         if (size_ > 0) {
            for (int i = 0; i < size_; ++i) {
               ptrs_[i] = other.ptrs_[i];
            }
         }
         // Nullify unused elements of ptrs_ array
         if (capacity_ > size_) {
            for (int i = size_; i < capacity_; ++i) {
               ptrs_[i] = 0;
            }
         }
      }
   }

   /*
   * Destructor.
   */
   template <typename Data>
   GStack<Data>::~GStack()
   {
      size_ = 0; 
      if (isAllocated()) {
         Memory::deallocate(ptrs_, capacity_);
         capacity_ = 0; 
      }
   }

   /*
   * Assignment, element by element.
   */
   template <typename Data>
   GStack<Data>& GStack<Data>::operator=(const GStack<Data>& other) 
   {
      // Check for self assignment
      if (this == &other) return *this;

      clear();
      for (int i = 0; i < other.size_; ++i) {
         append(other[i]);
      }
      return *this;
   }

   /*
   * Reserve space for the underlying array of Data* pointers.
   */
   template <typename Data>
   void GStack<Data>::reserve(int capacity) 
   {
      if (capacity <= 0) {
         UTIL_THROW("Cannot reserve with capacity <=0");
      }
      if (ptrs_ == 0) {
         assert(capacity_ == 0);
         assert(size_ == 0);
         Memory::allocate<Data*>(ptrs_, capacity);
         capacity_ = capacity;
         size_ = 0;
         for (int i = 0; i < capacity_; ++i) {
            ptrs_[i] = 0;
         }
      } else if (capacity > capacity_) {
         assert(capacity_ > 0);
         assert(capacity_ >= size_);
         assert(size_ >= 0);
         Data** newPtr = 0;
         Memory::allocate<Data*>(newPtr, capacity);
         if (size_ > 0) {
            for (int i = 0; i < size_; ++i) {
               newPtr[i] = ptrs_[i];
            }
         }
         if (size_ < capacity) {
            for (int i = size_; i < capacity; ++i) {
               newPtr[i] = 0;
            }
         }
         Memory::deallocate<Data*>(ptrs_, capacity_);
         ptrs_ = newPtr;
         capacity_ = capacity;
      }
   }

   /*
   * Deallocate associated memory.
   */
   template <typename Data>
   void GStack<Data>::deallocate() 
   {  
      size_ = 0; 
      if (isAllocated()) {
         assert(capacity_ > 0);
         Memory::deallocate<Data*>(ptrs_, capacity_);
         capacity_ = 0; 
      } 
   }

   /*
   * Reset to empty state, without deallocating.
   */
   template <typename Data>
   inline void GStack<Data>::clear()
   {  size_ = 0; }

   /*
   * Push an element onto the GStack.
   */
   template <typename Data>
   void GStack<Data>::push(Data& data) 
   {
      if (!isAllocated()) {
         assert(capacity_ == 0);
         Memory::allocate<Data*>(ptrs_, 64);
         capacity_ = 64;
         size_ = 0;
         for (int i = 0; i < capacity_; ++i) {
            ptrs_[i] = 0;
         }
      } else if (size_ == capacity_) {
         assert(capacity_ > 0);
         assert(capacity_ >= size_);
         assert(size_ >= 0);
         Data** newPtr = 0;
         Memory::allocate<Data*>(newPtr, 2*capacity_);
         if (size_ > 0) {
            for (int i = 0; i < size_; ++i) {
               newPtr[i] = ptrs_[i];
            }
         }
         if (size_ < 2*capacity_) {
            for (int i = size_; i < 2*capacity_; ++i) {
               newPtr[i] = 0;
            }
         }
         Memory::deallocate<Data*>(ptrs_, capacity_);
         ptrs_ = newPtr;
         capacity_ = 2*capacity_;
         // size_ is unchanged
      }

      // Append new element
      assert(size_ <= capacity_);
      ptrs_[size_] = &data;
      ++size_;
   }

   /*
   * Pop a Data object off the stack.
   */
   template <typename Data>
   Data& GStack<Data>::pop()
   {
      if (size_ == 0) {
         UTIL_THROW("Attempt to pop from empty stack");
      }
      Data *ptr = ptrs_[size_-1];
      ptrs_[size_-1] = 0;
      --size_;
      return *ptr;
   }

   /* 
   * Return a reference to the top element, without popping.
   */
   template <typename Data>
   inline Data& GStack<Data>::peek()
   {  return *ptrs_[size_-1]; }

   /* 
   * Return a const reference to the top element, without popping.
   */
   template <typename Data>
   inline const Data& GStack<Data>::peek() const
   {  return *ptrs_[size_-1]; }

   /*
   * Return current capacity.
   */
   template <typename Data>
   inline int GStack<Data>::capacity() const
   {  return capacity_; }

   /*
   * Return logical size.
   */
   template <typename Data>
   inline int GStack<Data>::size() const
   {  return size_; }

   /*
   * Is this GStack allocated?
   */
   template <class Data>
   inline bool GStack<Data>::isAllocated() const
   {  return (bool)ptrs_; }

   /* 
   * Return true if valid, or throw an exception.
   */
   template <typename Data>
   bool GStack<Data>::isValid() const
   {
      if (capacity_ < 0) {
         UTIL_THROW("Negative capacity_");
      }
      if (size_ < 0) {
         UTIL_THROW("Negative size_");
      }
      if (size_ > capacity_) {
         UTIL_THROW("size_ > capacity_");
      }

      if (ptrs_ != 0) {
         int i;
         for (i = 0; i < size_ ; ++i) {
            if (ptrs_[i] == 0) {
               UTIL_THROW("Null ptrs_[i] for i < size_");
            }
         }
         for (i = size_; i < capacity_ ; ++i) {
            if (ptrs_[i] != 0) {
               UTIL_THROW("Non-null ptrs_[i] for i >= size_");
            }
         }
      } else {
         if (capacity_ != 0) {
            UTIL_THROW("Unallocated stack but capacity != 0");
         }
      }
      return true;
   }

} 
#endif
