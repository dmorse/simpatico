#ifndef UTIL_CARRAY_2D_PARAM_H
#define UTIL_CARRAY_2D_PARAM_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Parameter.h>
#include <util/global.h>

#ifdef UTIL_MPI
#include <util/mpi/MpiSendRecv.h>
#endif

#include <iomanip> 

namespace Util
{

   /**
   * A Parameter associated with a 2D built-in C array.
   *
   * \ingroup Param_Module
   */
   template <class Type>
   class CArray2DParam : public Parameter
   {
      
   public:
  
      /** 
      * Constructor.
      *
      * Example: A 2 X 2 matrix stored in an oversized 3 x 3 C array.
      * \code
      *    double                matrix[3][3];
      *    CArray2DParam<double> param("matrix", matrix[0], 2, 2, 3);
      * \endcode
      *
      * \param label  parameter label (usually a literal C-string)
      * \param value  pointer to first element of associated 1D array
      * \param m      logical number of rows
      * \param n      logical number of columns
      * \param np     physical number of columns (elements per row).
      */
      CArray2DParam(const char *label, Type* value, int m, int n, int np);
  
      /** 
      * Read 2D array from file.
      */
      void readParam(std::istream &in);
  
      /**
      * Write 2D C array to file.
      */ 
      void writeParam(std::ostream &out);

      /**
      * Load 2D C array from archive.
      *
      * \param ar loading (input) archive.
      */
      void load(Serializable::IArchive& ar);

      /**
      * Save 2D C array to an archive.
      *
      * \param ar saving (output) archive.
      */
      void save(Serializable::OArchive& ar);

   protected:
   
      /// Pointer to first element in associated 1D C array
      Type* value_;
   
      /// Number of rows in array[][np]
      int m_; 

      /// Logical number of columns in array[][np]
      int n_; 

      /// Physical number of columns in array[][np]
      int np_; 
   
   
   };


   /*
   * CArray2D constructor.
   */
   template <class Type>
   CArray2DParam<Type>::CArray2DParam(const char* label, Type* value, int m, int n, int np)
    : Parameter(label),
      value_(value),
      m_(m),
      n_(n),
      np_(n)
   {}

   /*
   * Read CArray2DParam.
   */
   template <class Type>
   void CArray2DParam<Type>::readParam(std::istream &in)
   {
      if (isIoProcessor()) {
         int i, j;
         in >> label_;
         for (i = 0; i < m_; ++i) {
            for (j = 0; j < n_; ++j) {
               in >> value_[i*np_ + j];
            }
         }
         if (ParamComponent::echo()) {
            writeParam(Log::file());
         }
      }
      #ifdef UTIL_MPI
      if (hasIoCommunicator()) {
         bcast<Type>(ioCommunicator(), value_, m_*np_, 0); 
      }
      #endif
   }

   /*
   * Write a CArray2DParam.
   */
   template <class Type>
   void CArray2DParam<Type>::writeParam(std::ostream &out)
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
                << value_[i*np_ + j];
         }
         out << std::endl;
      }
   }

   /*
   * Load from an archive.
   */
   template <class Type>
   void CArray2DParam<Type>::load(Serializable::IArchive& ar)
   {
      if (isIoProcessor()) {
         int i, j;
         for (i = 0; i < m_; ++i) {
            for (j = 0; j < n_; ++j) {
               ar >> value_[i*np_ + j];
            }
         }
         if (ParamComponent::echo()) {
            writeParam(Log::file());
         }
      }
      #ifdef UTIL_MPI
      if (hasIoCommunicator()) {
         // Broadcast block containing m_ rows, each of np_ elements 
         bcast<Type>(ioCommunicator(), value_, m_*np_, 0); 
      }
      #endif
   }

   /*
   * Save to an archive.
   */
   template <class Type>
   void CArray2DParam<Type>::save(Serializable::OArchive& ar)
   {
      int i, j;
      for (i = 0; i < m_; ++i) {
         for (j = 0; j < n_; ++j) {
            ar << value_[i*np_ + j];
         }
      }
   }

} 
#endif
