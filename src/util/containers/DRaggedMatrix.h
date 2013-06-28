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

      /**
      * Destructor.
      *
      * Delete dynamically allocated C array.
      */
      ~DRaggedMatrix();

      /**
      * Allocate memory for a ragged matrix.
      */
      void allocate(const DArray<int>& rowSizes);

      /**
      * Return true iff this DRaggedMatrix has been allocated.
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

   /*
   * Return true if the DRaggedMatrix has been allocated, false otherwise.
   */
   template <class Data>
   inline bool DRaggedMatrix<Data>::isAllocated() const 
   {  return !(data_ == 0); }

}
#endif
