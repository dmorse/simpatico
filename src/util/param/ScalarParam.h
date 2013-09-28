#ifndef UTIL_SCALAR_PARAM_H
#define UTIL_SCALAR_PARAM_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Parameter.h>
#include <util/archives/Serializable_includes.h>
#include <util/global.h>

#ifdef UTIL_MPI
#include <util/mpi/MpiSendRecv.h>
#endif

#include <iomanip> 

namespace Util
{

   /** 
   * Template for a Parameter object associated with a scalar variable.
   *
   * This template can be used to define a Parameter subclass for any 
   * data type for which there exist inserter (<<) and extractor (>>) 
   * operators for stream io. 
   *
   * \ingroup Param_Module
   */
   template <class Type>
   class ScalarParam : public Parameter
   {
   
   public:
  
      /**
      * Constructor.
      *
      * \param label  label string const.
      * \param value  reference to parameter value.
      * \param isRequired  Is this a required parameter?
      */
      ScalarParam(const char *label, Type& value, bool isRequired = true); 
 
      /** 
      * Read parameter from stream.
      *
      * \param in input stream
      */
      void readParam(std::istream& in); 
 
      /**
      * Load from an archive.
      *
      * \param ar loading (input) archive.
      */
      void load(Serializable::IArchive& ar);

      /** 
      * Write parameter to stream.
      *
      * \param out output stream
      */
      void writeParam(std::ostream& out);

      /**
      * Save to an archive.
      *
      * \param ar saving (output) archive.
      */
      void save(Serializable::OArchive& ar);

      /**
      * Set the pointer to point a specific variable.
      *
      * \param value variable that holds the parameter value. 
      */
      void setValue(Type& value);

   protected:
   
      /// Pointer to value.
      Type* valuePtr_;

      /// Is value initialized?
      bool isInitialized_;

   private:

      /// Private and not implemented to prevent copying.
      ScalarParam(const ScalarParam<Type>& other);

      /// Private and not implemented to prevent assignment.
      ScalarParam<Type> operator = (const ScalarParam<Type>& other);
   
   };

   // ScalarParam Methods

   /*
   * ScalarParam<Type> constructor.
   */
   template <class Type>
   ScalarParam<Type>::ScalarParam(const char *label, Type &value, 
                                  bool isRequired)
    : Parameter(label, isRequired),
      valuePtr_(&value),
      isInitialized_(false)
   {}

   /*
   * Read a parameter.
   */
   template <class Type>
   void ScalarParam<Type>::readParam(std::istream &in)
   {
      if (isIoProcessor()) {
         in >> label_;
         if (Label::isClear()) {
            in >> *valuePtr_;
            isInitialized_ = true;
            if (ParamComponent::echo()) {
               writeParam(Log::file());
            }
         }
      }
      #ifdef UTIL_MPI
      if (hasIoCommunicator()) {
         is (isRequired()) {
            assert(success);
            bcast<Type>(ioCommunicator(), *valuePtr_, 0); 
         } else {
            bcast<bool>(ioCommunicator(), isInitialized_, 0); 
            if (success) {
               bcast<Type>(ioCommunicator(), *valuePtr_, 0); 
            }
         }
      }
      #endif
   }

   /*
   * Write a parameter.
   */
   template <class Type>
   void ScalarParam<Type>::writeParam(std::ostream& out)
   {
      if (isInitialized_) {
         out << indent();
         out << label_;
         out << std::right << std::scientific 
             << std::setprecision(Parameter::Precision) 
             << std::setw(Parameter::Width) << *valuePtr_ 
             << std::endl;
      }
   }

   /*
   * Load from an archive.
   */
   template <class Type>
   void ScalarParam<Type>::load(Serializable::IArchive& ar)
   {
      if (isIoProcessor()) {
         if (isRequired()) {
            isInitialized_ = true;
         } else {
            ar & isInitialized_;
         }
         if (isInitialized_) {
            if (valuePtr_) {  
               ar & *valuePtr_; 
            } else {
               UTIL_THROW("Attempt to load into null valuePtr_");
            }
            if (ParamComponent::echo()) {
               writeParam(Log::file());
            }
         }
      }
      #ifdef UTIL_MPI
      if (hasIoCommunicator()) {
         if (isRequired()) {
            isInitialized_ = true;
         } else {
            bcast<bool>(ioCommunicator(), isInitialized_, 0); 
         }
         if (isInitialized_) {
            bcast<Type>(ioCommunicator(), *valuePtr_, 0); 
         }
      }
      #endif
   }

   /*
   * Save to an archive.
   */
   template <class Type>
   void ScalarParam<Type>::save(Serializable::OArchive& ar)
   {
      if (!isRequired()) {
        ar & isInitialized_;
      }
      if (isInitialized_) {
         if (valuePtr_) {
            ar & *valuePtr_; 
         } else {
            UTIL_THROW("Attempt to write from null valuePtr_");
         }
      }
   }

   /*
   * Set the pointer to the parameter value.
   */
   template <class Type>
   void ScalarParam<Type>::setValue(Type& value)
   {  valuePtr_ = &value; }

} 
#endif
