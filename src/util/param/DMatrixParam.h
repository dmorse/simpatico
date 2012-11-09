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

class ParamTest;

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
      * Example: 
      * \code
      *    double                matrix[3][2];
      *    DMatrixParam<double> param("matrix", matrix[0], 3, 2);
      * \endcode
      *
      * \param label  parameter label (usually a literal C-string)
      * \param matrix DMatrix object
      * \param m      number of rows
      * \param n      number of columns
      */
      DMatrixParam(const char *label, DMatrix<Type>& matrix, int m, int n);
  
      /** 
      * Read 2D array from file.
      */
      void readParam(std::istream &in);
  
      /**
      * Write 2D C array to file.
      */ 
      void writeParam(std::ostream &out);

   protected:
   
      /// Pointer to first element in associated 1D C array
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
      if (isParamIoProcessor()) {
         int i, j;
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
      if (hasParamCommunicator()) {
         bcast<Type>(paramCommunicator(), *matrixPtr_, m_, n_, 0); 
      }
      #endif
   }

   /*
   * Write a DMatrixParam.
   */
   template <class Type>
   void DMatrixParam<Type>::writeParam(std::ostream &out)
   {
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

} 
#endif
