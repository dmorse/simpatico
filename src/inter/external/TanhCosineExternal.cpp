#ifndef INTER_TANH_COSINE_EXTERNAL_CPP
#define INTER_TANH_COSINE_EXTERNAL_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, Jian Qin and David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "TanhCosineExternal.h"

#include <iostream>

namespace Inter
{

   using namespace Util;

   /* 
   * Constructor.
   */
   TanhCosineExternal::TanhCosineExternal() 
    : perpDirection_(0),
      width_(),
      externalParameter_(),
      periodicity_(),
      boundaryPtr_(0),
      nAtomType_(0), 
      isInitialized_(false)
   { setClassName("TanhCosineExternal"); }
   
   /* 
   * Copy constructor.
   */
   TanhCosineExternal::TanhCosineExternal(const TanhCosineExternal& other)
    : perpDirection_(other.perpDirection_),
      width_(other.width_),
      externalParameter_(other.externalParameter_),
      periodicity_(other.periodicity_),
      boundaryPtr_(other.boundaryPtr_),
      nAtomType_(other.nAtomType_),
      isInitialized_(other.isInitialized_)
   {
      prefactor_.allocate(nAtomType_);
      for (int i=0; i < nAtomType_; ++i) {
        prefactor_[i] = other.prefactor_[i];
      }
   } 
     
   /* 
   * Assignment operator.
   */
   TanhCosineExternal& TanhCosineExternal::operator = (const TanhCosineExternal& other)
   {
      perpDirection_    = other.perpDirection_;
      boundaryPtr_      = other.boundaryPtr_;
      nAtomType_        = other.nAtomType_;
      isInitialized_    = other.isInitialized_;
      width_            = other.width_;
      periodicity_      = other.periodicity_;
      externalParameter_   = other.externalParameter_;
      int i;
      for (i=0; i < nAtomType_; ++i) {
        prefactor_[i] = other.prefactor_[i];
      }
      return *this;
   }

   /* 
   * Set nAtomType
   */
   void TanhCosineExternal::setNAtomType(int nAtomType) 
   {  
      if (nAtomType <= 0) {
         UTIL_THROW("nAtomType <= 0");
      }
      if (nAtomType > MaxAtomType) {
         UTIL_THROW("nAtomType > TanhCosineExternal::MaxAtomType");
      }
      nAtomType_ = nAtomType;
   }

   void TanhCosineExternal::setExternalParameter(double externalParameter) 
   {  
      // Preconditions
      if (!isInitialized_) {
         UTIL_THROW("TanhCosineExternal potential is not initialized");
      }
     
      externalParameter_ = externalParameter;
   }

   
   /* 
   * Set pointer to the Boundary.
   */
   void TanhCosineExternal::setBoundary(Boundary &boundary)
   {  boundaryPtr_ = &boundary; }
   
   /* 
   * Read potential parameters from file.
   */
   void TanhCosineExternal::readParameters(std::istream &in) 
   {
      if (nAtomType_ == 0) {
         UTIL_THROW("nAtomType must be set before readParam");
      }
      if (!boundaryPtr_) {
         UTIL_THROW("Boundary must be set before readParam");
      }
   
      // Read parameters
      read<int>(in, "perpDirection", perpDirection_);
      if (perpDirection_ < 0 || perpDirection_ >= Dimension) {
         UTIL_THROW("Invalid index for perpendicular direction.");
      }
      prefactor_.allocate(nAtomType_);
      readDArray<double>(in, "prefactor", prefactor_, nAtomType_);
      read<double>(in, "externalParameter", externalParameter_);
      read<double>(in, "interfaceWidth", width_);
      read<int>(in, "periodicity", periodicity_);

      isInitialized_ = true;
   }

   /*
   * Load internal state from an archive.
   */
   void TanhCosineExternal::loadParameters(Serializable::IArchive &ar)
   {
      ar >> nAtomType_; 
      if (nAtomType_ <= 0) {
         UTIL_THROW( "nAtomType must be positive");
      }
      loadParameter<int>(ar, "perpDirection", perpDirection_);
      if (perpDirection_ < 0 || perpDirection_ >= Dimension) {
         UTIL_THROW("Invalid index for perpendicular direction.");
      }
      prefactor_.allocate(nAtomType_);
      loadDArray<double>(ar, "prefactor", prefactor_, nAtomType_);
      loadParameter<double>(ar, "externalParameter", externalParameter_);
      loadParameter<double>(ar, "interfaceWidth", width_);
      loadParameter<int>(ar, "periodicity", periodicity_);
      isInitialized_ = true;
   }

   /*
   * Save internal state to an archive.
   */
   void TanhCosineExternal::save(Serializable::OArchive &ar)
   {
      ar << nAtomType_;
      ar << perpDirection_;
      ar << externalParameter_;
      ar << prefactor_;
      ar << width_;
      ar << periodicity_;
   }

   double TanhCosineExternal::externalParameter() const
   {  return externalParameter_; }

   /*
   * Return name string "TanhCosineExternal".
   */
   std::string TanhCosineExternal::className() const
   {  return std::string("TanhCosineExternal"); }
 
} 
#endif
