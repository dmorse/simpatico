#ifndef UTIL_D_RAGGED_MATRIX_H
#define UTIL_D_RAGGED_MATRIX_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/RaggedMatrix.h>
#include <util/containers/DArray.h>
#include <util/global.h>

namespace Util
{

   /**
   * Dynamically allocated RaggedMatrix.
   *
   * \ingroup Matrix_Module
   */
   template <typename Data>
   class DRaggedMatrix : public RaggedMatrix<Data>
   {

      using RaggedMatrix<Data>::data_;
      using RaggedMatrix<Data>::rows_;
      using RaggedMatrix<Data>::capacity1_;
      using RaggedMatrix<Data>::capacity2_;

   public:

      /**
      * Constructor.
      */
      DRaggedMatrix();

      #if 0
      /**
      * Copy constructor.
      */
      DRaggedMatrix(const DRaggedMatrix<Data>& other);
      #endif

      /**
      * Destructor.
      *
      * Delete dynamically allocated C array.
      */
      ~DRaggedMatrix();

      #if 0
      /**
      * Assignment.
      *
      * \throw Exception if LHS and RHS dimensions do not match.
      */
      DRaggedMatrix<Data>& operator= (const DRaggedMatrix<Data>& other);
      #endif

      /**
      * Allocate memory for a ragged matrix.
      */
      void allocate(const DArray<int>& rowSizes);

      /**
      * Return true if the DRaggedMatrix has been allocated, false otherwise.
      */
      bool isAllocated() const;
   

   };

   // Method definitions

   /*
   * Default constructor.
   */
   template <typename Data>
   DRaggedMatrix<Data>::DRaggedMatrix() :
      RaggedMatrix<Data>()
   {}
     
   #if 0 
   /*
   * Copy constructor.
   */
   template <typename Data>
   DRaggedMatrix<Data>::DRaggedMatrix(const DRaggedMatrix<Data>& other) 
     : RaggedMatrix<Data>()
   {

      // Precondition
      if (other.data_ == 0) {
         UTIL_THROW("Other DRaggedMatrix must be allocated");
      }

      allocate(other.capacity1_, other.capacity2_);

      // Copy elements
      for (int i = 0; i < capacity1_*capacity2_; ++i) {
         data_[i] = other.data_[i];
      }
   }
   #endif
  
   #if 0 
   /*
   * Assignment.
   */
   template <typename Data>
   DRaggedMatrix<Data>& DRaggedMatrix<Data>::operator = (const DRaggedMatrix<Data>& other) 
   {

      // Check for self assignment.
      if (this == &other) return *this;

      // Precondition
      if (other.data_ == 0) {
         UTIL_THROW("RHS DRaggedMatrix must be allocated in assignment");
      }

      // If this this DRaggedMatrix if not allocated, allocate now.
      // If this this DRaggedMatrix is allocated, check that capacities are equal.
      if (data_ == 0) {
         allocate(other.capacity1_, other.capacity2_);
      } else {
         if (capacity1_ != other.capacity1_ || 
             capacity2_ != other.capacity2_) {
            UTIL_THROW("Unequal capacities in assignment");
         }
      }

      // Copy elements
      for (int i = 0; i < capacity1_*capacity2_; ++i) {
         data_[i] = other.data_[i];
      }

      return *this;
   }
   #endif

   /*
   * Destructor.
   */
   template <typename Data>
   DRaggedMatrix<Data>::~DRaggedMatrix()
   {
      if (data_ != 0) {
         delete [] data_;
      }
      if (rows_ != 0) {
         delete [] rows_;
      }
      if (capacity2_ != 0) {
         delete [] capacity2_;
      }
   }

   /*
   * Allocate memory and set row sizes.
   */
   template <typename Data>
   void DRaggedMatrix<Data>::allocate(const DArray<int>& rowSizes)
   {
      int i, capacity;

      if (data_ != 0) 
         UTIL_THROW("Attempt to re-allocate a RaggedMatrix");

      if (rowSizes.capacity() <= 0) 
         UTIL_THROW("rowSizes.capacity() must be positive");

      capacity = 0; 
      for (i = 0; i < rowSizes.capacity(); ++i) {
         if (rowSizes[i] < 0) 
            UTIL_THROW("rowSizes must all be nonnegative");
         capacity += rowSizes[i];
      }

      if (capacity == 0) 
         UTIL_THROW("Sum of row sizes must be positive");

      // Allocate all memory
      capacity1_  = rowSizes.capacity();
      capacity2_  = new int[capacity1_];
      rows_       = new Data*[capacity1_];
      data_       = new Data[capacity];

      // Set row sizes and pointers to rows
      Data* ptr = data_;
      for (i = 0; i < capacity1_; ++i) {
         capacity2_[i] = rowSizes[i];
         rows_[i] = ptr;
         ptr += rowSizes[i];
      }
      
   }

   /// Return true if the DRaggedMatrix has been allocated, false otherwise.
   template <class Data>
   inline bool DRaggedMatrix<Data>::isAllocated() const 
   {  return !(data_ == 0); }

}
#endif
