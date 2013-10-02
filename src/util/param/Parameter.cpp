#ifndef UTIL_PARAMETER_CPP
#define UTIL_PARAMETER_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Parameter.h"

#ifdef UTIL_MPI
#include <util/mpi/MpiSendRecv.h>
#endif

namespace Util
{

   /*
   * Constructor. 
   */
   Parameter::Parameter(const char *label, bool isRequired)
    : label_(label, isRequired),
      isActive_(isRequired)
   {}

   /*
   * Destructor.
   */
   Parameter::~Parameter()
   {}

   /*
   * Read a parameter.
   */
   void Parameter::readParam(std::istream &in)
   {
      if (isIoProcessor()) {
         in >> label_;
         if (Label::isClear()) {
            readValue(in);
            isActive_ = true;
            if (ParamComponent::echo()) {
               writeParam(Log::file());
            }
         } else {
            if (ParamComponent::echo()) {
               Log::file() << indent() 
                           << "[ " << label_.string() 
                           << " is absent ]" << std::endl;
            }
         }
      }
      #ifdef UTIL_MPI
      if (hasIoCommunicator()) {
         if (isRequired()) {
            bcastValue();
            isActive_ = true;
         } else {
            bcast<bool>(ioCommunicator(), isActive_, 0); 
            if (isActive_) {
               bcastValue();
            }
         }
      }
      #endif
   }

   /*
   * Load from an archive.
   */
   void Parameter::load(Serializable::IArchive& ar)
   {
      if (isIoProcessor()) {
         if (!isRequired()) {
            ar & isActive_;
         }
         if (isActive_) {
            loadValue(ar);
            if (ParamComponent::echo()) {
               writeParam(Log::file());
            }
         }
      }
      #ifdef UTIL_MPI
      if (hasIoCommunicator()) {
         if (!isRequired()) {
            bcast<bool>(ioCommunicator(), isActive_, 0); 
         }
         if (isActive_) {
            bcastValue();
         }
      }
      #endif
   }

   /*
   * Save to an archive.
   */
   void Parameter::save(Serializable::OArchive& ar)
   {
      if (!isRequired()) {
        ar & isActive_;
      }
      if (isActive_) {
         saveValue(ar);
      }
   }

   /*
   * Return label string.
   */
   std::string Parameter::label() const
   {  return label_.string(); }

   /*
   * Is this a required parameter?
   */
   bool Parameter::isRequired() const
   {  return label_.isRequired(); }

   /*
   * Is this an active parameter?
   */
   bool Parameter::isActive() const
   {  return isActive_; }

} 
#endif
