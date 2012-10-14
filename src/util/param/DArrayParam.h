#ifndef UTIL_D_ARRAY_PARAM_H
#define UTIL_D_ARRAY_PARAM_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Parameter.h>    // base class
#include <util/containers/DArray.h>  // member
#include <util/global.h>

#include <iomanip> 

class ParamTest;

namespace Util
{

   /**
   * A Parameter associated with a DArray container.
   *
   * \ingroup Param_Module
   */
   template <class Type>
   class DArrayParam : public Parameter
   {
   
   public:
   
      /// Constructor.
      DArrayParam(const char *label, DArray<Type>& array, int n = 0);
 
      /** 
      * Read parameter from stream.
      *
      * \param in input stream
      */
      void readParam(std::istream &in);
 
      /** 
      * Write parameter to stream.
      *
      * \param out output stream
      */
      void writeParam(std::ostream &out);

   protected:
   
      /// Pointer to associated DArray.
      DArray<Type>* arrayPtr_;
   
      /// Logical array dimension
      int n_;
   
   };

   /*
   * DArrayParam<Type> constructor.
   */
   template <class Type>
   DArrayParam<Type>::DArrayParam(const char *label, DArray<Type>& array, int n)
    : Parameter(label),
      arrayPtr_(&array),
      n_(n)
   {}

   /*
   * Read a DArray parameter.
   */
   template <class Type>
   void DArrayParam<Type>::readParam(std::istream &in)
   {

      if (!(arrayPtr_->isAllocated())) {
         UTIL_THROW("Cannot read unallocated DArray");
      }
      if (arrayPtr_->capacity() < n_) {
         UTIL_THROW("Error: DArray capacity < n");
      }

      if (isParamIoProcessor()) {
         in >> label_;
         for (int i = 0; i < n_; ++i) {
            in >> (*arrayPtr_)[i];
         }
         if (ParamComponent::echo()) {
            writeParam(Log::file());
         }
      }
      #ifdef UTIL_MPI
         if (hasParamCommunicator()) 
            bcast<Type>(paramCommunicator(), *arrayPtr_, n_, 0); 
      #endif

   }

   /*
   * Write a DArray parameter.
   */
   template <class Type>
   void DArrayParam<Type>::writeParam(std::ostream &out) 
   {

      // Preconditions
      if (!(arrayPtr_->isAllocated())) {
         UTIL_THROW("Cannot read unallocated DArray");
      }
      if (arrayPtr_->capacity() < n_) {
         UTIL_THROW("Error: DArray capacity < n");
      }

      Label space("");
      int i;

      for (i = 0; i < n_; ++i) {
         if (i == 0) {
            out << indent() << label_;
         } else {
            out << indent() << space;
         }
         out << std::right << std::scientific 
             << std::setprecision(Parameter::Precision) 
             << std::setw(Parameter::Width)
             << (*arrayPtr_)[i] 
             << std::endl;
      }

   }

} 
#endif
