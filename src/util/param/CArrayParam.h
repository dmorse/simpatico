#ifndef UTIL_CARRAY_PARAM_H
#define UTIL_CARRAY_PARAM_H

#include <util/param/Parameter.h>
#include <util/global.h>

#ifdef UTIL_MPI
#include <util/mpi/MpiSendRecv.h>
#endif

#include <iomanip> 

namespace Util
{

   /**
   * A Parameter associated with a 1D C array. 
   *
   * \ingroup Param_Module
   */
   template <class Type>
   class CArrayParam : public Parameter
   {
   
   public:
   
      /// Constructor.
      CArrayParam(const char *label, Type *value, int n=0);
 
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

      /**
      * Load from an archive.
      *
      * \param ar loading (input) archive.
      */
      void load(Serializable::IArchive& ar);

      /**
      * Save to an archive.
      *
      * \param ar saving (output) archive.
      */
      void save(Serializable::OArchive& ar);

 
   protected:
   
      /// Pointer to value.
      Type* value_;
   
      /// Array dimension
      int n_;    
   
   };

   /*
   * CArrayParam<Type> constructor.
   */
   template <class Type>
   CArrayParam<Type>::CArrayParam(const char *label, Type* value, int n)
    : Parameter(label),
      value_(value),
      n_(n)
   {}

   /*
   * Read a C array parameter.
   */
   template <class Type>
   void CArrayParam<Type>::readParam(std::istream &in)
   {
      if (isIoProcessor()) {
         in >> label_;
         in >> value_[0];
         for (int i = 1; i < n_; ++i) {
            in >> value_[i];
         }
         if (ParamComponent::echo()) {
            writeParam(Log::file());
         }
      }
      #ifdef UTIL_MPI
      if (hasIoCommunicator()) {
         bcast<Type>(ioCommunicator(), value_, n_, 0); 
      }
      #endif
   }

   /*
   * Write a C array
   */
   template <class Type>
   void CArrayParam<Type>::writeParam(std::ostream &out) 
   {

      #if 0
      // Output label
      out << indent();
      out << label_;

      // Output values
      for (int i = 0; i < n_; ++i) {
         out.setf(std::ios::scientific);
         out.width(Parameter::Width);
         out.precision(Parameter::Precision);
         out << value_[i];
      }
      out << std::endl;
      #endif

      Label space("");
      for (int i = 0; i < n_; ++i) {
         if (i == 0) {
            out << indent() << label_;
         } else {
            out << indent() << space;
         }
         out << std::right << std::scientific 
             << std::setprecision(Parameter::Precision) 
             << std::setw(Parameter::Width)
             << value_[i] 
             << std::endl;
      }

   }

   /*
   * Load from an archive.
   */
   template <class Type>
   void CArrayParam<Type>::load(Serializable::IArchive& ar)
   {
      if (isIoProcessor()) {
         for (int i = 0; i < n_; ++i) {
            ar >> value_[i];
         }
         if (ParamComponent::echo()) {
            writeParam(Log::file());
         }
      }
      #ifdef UTIL_MPI
      if (hasIoCommunicator()) {
         bcast<Type>(ioCommunicator(), value_, n_, 0); 
      }
      #endif
   }

   /*
   * Save to an archive.
   */
   template <class Type>
   void CArrayParam<Type>::save(Serializable::OArchive& ar)
   {
      for (int i = 0; i < n_; ++i) {
         ar << value_[i];
      }
   }

} 
#endif
