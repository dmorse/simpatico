/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "LocalLamellarOrderingExternal.h"

#include <iostream>

namespace Simp
{

   using namespace Util;

   /* 
   * Constructor.
   */
   LocalLamellarOrderingExternal::LocalLamellarOrderingExternal() 
    : perpDirection_(0),
      parallelDirection_(0),
      fraction_(0.0),
      width_(),
      externalParameter_(),
      periodicity_(),
      boundaryPtr_(0),
      nAtomType_(0), 
      isInitialized_(false)
   { setClassName("LocalLamellarOrderingExternal"); }
   
   /* 
   * Copy constructor.
   */
   LocalLamellarOrderingExternal::LocalLamellarOrderingExternal(const LocalLamellarOrderingExternal& other)
    : perpDirection_(other.perpDirection_),
      parallelDirection_(other.parallelDirection_),
      fraction_(other.fraction_),
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
   LocalLamellarOrderingExternal& LocalLamellarOrderingExternal::operator = (const LocalLamellarOrderingExternal& other)
   {
      perpDirection_    = other.perpDirection_;
      parallelDirection_    = other.parallelDirection_;
      fraction_    = other.fraction_;
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
   void LocalLamellarOrderingExternal::setNAtomType(int nAtomType) 
   {  
      if (nAtomType <= 0) {
         UTIL_THROW("nAtomType <= 0");
      }
      if (nAtomType > MaxAtomType) {
         UTIL_THROW("nAtomType > LocalLamellarOrderingExternal::MaxAtomType");
      }
      nAtomType_ = nAtomType;
   }

   void LocalLamellarOrderingExternal::setExternalParameter(double externalParameter) 
   {  
      // Preconditions
      if (!isInitialized_) {
         UTIL_THROW("LocalLamellarOrderingExternal potential is not initialized");
      }
     
      externalParameter_ = externalParameter;
   }

   
   /* 
   * Set pointer to the Boundary.
   */
   void LocalLamellarOrderingExternal::setBoundary(Boundary &boundary)
   {  boundaryPtr_ = &boundary; }
   
   /* 
   * Read potential parameters from file.
   */
   void LocalLamellarOrderingExternal::readParameters(std::istream &in) 
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
      read<int>(in, "parallelDirection", parallelDirection_);
      if (parallelDirection_ < 0 || parallelDirection_ >= Dimension) {
         UTIL_THROW("Invalid index for parallel direction.");
      }
      read<double>(in, "fraction", fraction_);
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
   void LocalLamellarOrderingExternal::loadParameters(Serializable::IArchive &ar)
   {
      ar >> nAtomType_; 
      if (nAtomType_ <= 0) {
         UTIL_THROW( "nAtomType must be positive");
      }
      loadParameter<int>(ar, "perpDirection", perpDirection_);
      if (perpDirection_ < 0 || perpDirection_ >= Dimension) {
         UTIL_THROW("Invalid index for perpendicular direction.");
      }
      loadParameter<int>(ar, "parallelDirection", perpDirection_);
      if (parallelDirection_ < 0 || parallelDirection_ >= Dimension) {
         UTIL_THROW("Invalid index for parallel direction.");
      }
      loadParameter<double>(ar, "fraction", fraction_);
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
   void LocalLamellarOrderingExternal::save(Serializable::OArchive &ar)
   {
      ar << nAtomType_;
      ar << perpDirection_;
      ar << parallelDirection_;
      ar << fraction_;
      ar << externalParameter_;
      ar << prefactor_;
      ar << width_;
      ar << periodicity_;
   }

   /*
   * Set a potential energy parameter, identified by a string.
   */
   void LocalLamellarOrderingExternal::set(std::string name, double value)
   {
      if (name == "fraction") {
         fraction_ = value;
      } else
      if (name == "externalParameter") {
         externalParameter_ = value;
      } else
      if (name == "interfaceWidth") {
         width_ = value;
      } else {
         UTIL_THROW("Unrecognized parameter name");
      }
   }

   /*
   * Get a parameter value, identified by a string.
   */
   double LocalLamellarOrderingExternal::get(std::string name) const
   {
      double value = 0.0;
      if (name == "fraction") {
         value = fraction_;
      } else
      if (name == "externalParameter") {
         value = externalParameter_;
      } else
      if (name == "interfaceWidth") {
         value = width_;
      } else {
         UTIL_THROW("Unrecognized parameter name");
      }
      return value;
   }

   double LocalLamellarOrderingExternal::externalParameter() const
   {  return externalParameter_; }

   /*
   * Return name string "LocalLamellarOrderingExternal".
   */
   std::string LocalLamellarOrderingExternal::className() const
   {  return std::string("LocalLamellarOrderingExternal"); }
 
} 
