#ifndef UTIL_DMATRIX_PARAM_H
#define UTIL_DMATRIX_PARAM_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Parameter.h>
#include <util/containers/DMatrix.h>
#ifdef UTIL_MPI
#include <util/mpi/MpiSendRecv.h>
#endif
#include <util/global.h>

#include <iomanip> 

namespace Util
{

   /**
   * A Parameter associated with a 2D built-in C array.
   *
   * \ingroup Param_Module
   */
   template <class Type>
   class DMatrixParam : public Parameter
   {
      
   public:
  
      /** 
      * Constructor.
      *
      * \param label  parameter label (usually a literal C-string)
      * \param matrix DMatrix object
      * \param m      number of rows
      * \param n      number of columns
      */
      DMatrixParam(const char *label, DMatrix<Type>& matrix, int m, int n);
  
      /** 
      * Read DMatrix from file.
      */
      void readParam(std::istream &in);
  
      /**
      * Write DMatrix to file.
      */ 
      void writeParam(std::ostream &out);

      /**
      * Load DMatrix from archive.
      *
      * \param ar loading (input) archive.
      */
      void load(Serializable::IArchive& ar);

      /**
      * Save DMatrix to an archive.
      *
      * \param ar saving (output) archive.
      */
      void save(Serializable::OArchive& ar);


   protected:
   
      /// Pointer to associated DMatrix.
      DMatrix<Type>* matrixPtr_;
   
      /// Number of rows in array[m][n]
      int m_; 

      /// Number of columns in array[m][n]
      int n_; 
   
   };

   /*
   * DMatrix constructor.
   */
   template <class Type>
   DMatrixParam<Type>::DMatrixParam(const char* label, DMatrix<Type>& matrix, int m, int n)
    : Parameter(label),
      matrixPtr_(&matrix),
      m_(m),
      n_(n)
   {}

   /*
   * Read DMatrixParam.
   */
   template <class Type>
   void DMatrixParam<Type>::readParam(std::istream &in)
   {
      // Preconditions
      if (!(matrixPtr_->isAllocated())) {
         UTIL_THROW("Cannot read unallocated DMatrix");
      }
      if (m_ > matrixPtr_->capacity1()) {
         UTIL_THROW("Error: Logical size m_ > DMatrix<Type>::capacity1()");
      }
      if (n_ > matrixPtr_->capacity2()) {
         UTIL_THROW("Error: Logical size n_ > DMatrix<Type>::capacity2()");
      }

      // Read from file
      if (isIoProcessor()) {
         int i, j;
         //int m = matrixPtr_->capacity1();
         //int n = matrixPtr_->capacity2();
         in >> label_;
         for (i = 0; i < m_; ++i) {
            for (j = 0; j < n_; ++j) {
               in >> (*matrixPtr_)(i, j);
            }
         }
         if (ParamComponent::echo()) {
            writeParam(Log::file());
         }
      }

      #ifdef UTIL_MPI
      // Broadcast
      if (hasIoCommunicator()) {
         //int m = matrixPtr_->capacity1();
         //int n = matrixPtr_->capacity2();
         bcast<Type>(ioCommunicator(), *matrixPtr_, m_, n_, 0); 
      }
      #endif
   }

   /*
   * Write a DMatrixParam.
   */
   template <class Type>
   void DMatrixParam<Type>::writeParam(std::ostream &out)
   {
      // Preconditions
      if (!(matrixPtr_->isAllocated())) {
         UTIL_THROW("Cannot read unallocated DMatrix");
      }
      if (m_ > matrixPtr_->capacity1()) {
         UTIL_THROW("Error: Logical size m_ > DMatrix<Type>::capacity1()");
      }
      if (n_ > matrixPtr_->capacity2()) {
         UTIL_THROW("Error: Logical size n_ > DMatrix<Type>::capacity2()");
      }

      Label space("");
      int i, j;
      for (i = 0; i < m_; ++i) {
         if (i == 0) {
            out << indent() << label_;
         } else {
            out << indent() << space;
         }
         for (j = 0; j < n_; ++j) {
            out << std::right << std::scientific 
                << std::setprecision(Parameter::Precision) 
                << std::setw(Parameter::Width)
                << (*matrixPtr_)(i, j);
         }
         out << std::endl;
      }
   }

   /*
   * Load from an archive.
   */
   template <class Type>
   void DMatrixParam<Type>::load(Serializable::IArchive& ar)
   {
      if (!(matrixPtr_->isAllocated())) {
         matrixPtr_->allocate(m_, n_);
      }

      // Load from archive on ioProcessor
      if (isIoProcessor()) {
         ar >> *matrixPtr_;
         if (ParamComponent::echo()) {
            writeParam(Log::file());
         }
      }

      #ifdef UTIL_MPI
      // Broadcast to all processors
      if (hasIoCommunicator()) {
         int m = matrixPtr_->capacity1();
         int n = matrixPtr_->capacity2();
         bcast<Type>(ioCommunicator(), *matrixPtr_, m, n, 0); 
      }
      #endif
   }

   /*
   * Save to an archive.
   */
   template <class Type>
   void DMatrixParam<Type>::save(Serializable::OArchive& ar)
   {  ar << *matrixPtr_; }

} 
#endif
