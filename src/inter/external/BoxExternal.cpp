#ifndef BOX_EXTERNAL_CPP
#define BOX_EXTERNAL_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, Jian Qin and David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "BoxExternal.h"

#include <iostream>

namespace Inter
{

   using namespace Util;

   /* 
   * Constructor.
   */
   BoxExternal::BoxExternal() 
    : boundaryPtr_(0),
      nAtomType_(0),
      isInitialized_(false)
   {}
   
   /* 
   * Copy constructor.
   */
   BoxExternal::BoxExternal(const BoxExternal& other)
    : epsilon_(other.epsilon_),
      sigma_(other.sigma_),
      sigmaCb_(other.sigmaCb_),
      cutoff_(other.cutoff_),
      coeff_(other.coeff_),
      boundaryPtr_(other.boundaryPtr_),
      nAtomType_(other.nAtomType_),
      isInitialized_(other.isInitialized_)
   {}
   
   /* 
   * Assignment operator.
   */
   BoxExternal& BoxExternal::operator = (const BoxExternal& other)
   {
      boundaryPtr_   = other.boundaryPtr_;
      nAtomType_     = other.nAtomType_;
      isInitialized_ = other.isInitialized_;
      epsilon_       = other.epsilon_;
      sigma_         = other.sigma_;
      sigmaCb_       = other.sigmaCb_;
      cutoff_        = other.cutoff_;
      coeff_         = other.coeff_;
      return *this;
   }

   /* 
   * Set nAtomType
   */
   void BoxExternal::setNAtomType(int nAtomType) 
   {  
      if (nAtomType <= 0) {
         UTIL_THROW("nAtomType <= 0");
      }
      if (nAtomType > MaxAtomType) {
         UTIL_THROW("nAtomType > BoxExternal::MaxAtomType");
      }
      nAtomType_ = nAtomType;
   }

   /* 
   * Set pointer to the Boundary.
   */
   void BoxExternal::setBoundary(Boundary &boundary)
   {  boundaryPtr_ = &boundary; }

   /**
   * Set epsilon to the externalParameter.
   */
   void BoxExternal::setExternalParameter(double externalParameter)
   {
      // Preconditions
      if (!isInitialized_) {
         UTIL_THROW("BoxExternal potential is not initialized");
      }

      epsilon_ = externalParameter;
   }
   
   /* 
   * Read potential parameters from file.
   */
   void BoxExternal::readParameters(std::istream &in) 
   {
      if (nAtomType_ == 0) {
         UTIL_THROW("nAtomType must be set before readParam");
      }
      if (!boundaryPtr_) {
         UTIL_THROW("Boundary must be set before readParam");
      }
   
      // Read parameters
      // //readBegin(in,  "BoxExternal");
      read<double>(in, "epsilon", epsilon_);
      read<double>(in, "sigma", sigma_);
      read<double>(in, "cutoff", cutoff_);
      // //readEnd(in);

      // Initialize dependent variables (particle number density = 0.7)
      sigmaCb_ = sigma_ * sigma_ * sigma_;
      coeff_   = 4.0 * epsilon_ * acos(-1.0) / 45.0 * 0.7 * sigmaCb_;

      isInitialized_ = true;
   }

   /**
   * Return epsilon for externalParameter.
   */
   double BoxExternal::externalParameter() const
   {
     return epsilon_;
   }

   /**
   * Return name string "BoxExternal" for this evaluator class.
   */
   std::string BoxExternal::className() const
   {  return std::string("BoxExternal"); }
 

} 
#endif
