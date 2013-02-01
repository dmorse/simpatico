#ifndef UTIL_MATRIX_H
#define UTIL_MATRIX_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>

namespace Util
{

   /**
   * Two-dimensional array container template (abstract).
   *
   * An Matrix object A is a two-dimensional array in which the operator
   * A(i,j) returns a reference to element j of row i.
   *
   * The memory for a Matrix is stored in a single one-dimensional C array.
   *
   * A Matrix is abstract because it cannot allocate memory.
   * Concrete subclasses: DMatrix and FMatrix.
   *
   * \ingroup Matrix_Module
   */
   template <typename Data>
   class Matrix
   {

   public:

      // Protected default constructor, to prohibit direct instantiation.

      // Private copy constructor, to prohibit copy construction.

      /**
      * Destructor.
      */
      virtual ~Matrix();

      /**
      * Get number of rows.
      *
      * \return Number of columns (i.e., range of first array index)
      */
      int capacity1();

      /**
      * Get number of columns.
      *
      * \return Number of columns (i.e., range of second array index)
      */
      int capacity2();

      /**
      * Return element (i,j) of matrix by const reference.
      *
      * \param  i  row index.
      * \param  j  column index.
      */
      const Data& operator() (int i, int j) const;

      /**
      * Return element (i,j) of matrix by reference.
      *
      * \param  i  row index.
      * \param  j  column index.
      */
      Data& operator() (int i, int j);

   protected:

      /// Pointer to 1D C array of all elements.
      Data*  data_;

      /// Number of rows (range of first index).
      int  capacity1_;

      /// Number of columns (range of first index).
      int  capacity2_;

      /**
      * Default constructor. 
      *
      * Protected to prevent direct instantiation.
      */
      Matrix();

   private:

      /**
      * Copy constructor, private and not implemented.
      */
      Matrix(const Matrix<Data>& other);

      /**
      * Assignment, private and not implemented.
      */
      Matrix<Data>& operator = (const Matrix<Data>& other);

   }; 

   // Method definitions

   /**
   * Constructor (protected).
   */
   template <typename Data>
   inline Matrix<Data>::Matrix()
    : data_(0),
      capacity1_(0),
      capacity2_(0)
   {}

   /*
   * Destructor.
   */
   template <typename Data>
   inline Matrix<Data>::~Matrix()
   {}

   /*
   * Get number of rows.
   */
   template <typename Data>
   inline int Matrix<Data>::capacity1()
   { return capacity1_; }

   /*
   * Get number of columns.
   */
   template <typename Data>
   inline int Matrix<Data>::capacity2()
   { return capacity2_; }

   /*
   * Return element (i,j) of matrix by const reference.
   */
   template <typename Data>
   inline const Data& Matrix<Data>::operator() (int i, int j) const
   {
      assert(data_ != 0);
      assert(i >= 0);
      assert(i < capacity1_);
      assert(j >= 0);
      assert(j < capacity2_);
      return *(data_ + i*capacity2_ + j);
   }

   /*
   * Return element (i,j) of matrix by reference.
   */
   template <typename Data>
   inline Data& Matrix<Data>::operator() (int i, int j)
   {
      assert(data_ != 0);
      assert(i >= 0);
      assert(i < capacity1_);
      assert(j >= 0);
      assert(j < capacity2_);
      return *(data_ + i*capacity2_ + j);
   }


}
#endif
