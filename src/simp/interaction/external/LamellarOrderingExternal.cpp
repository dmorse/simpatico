/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "LamellarOrderingExternal.h"

#include <iostream>

namespace Simp
{

   using namespace Util;

   /*
   * Constructor.
   */
   LamellarOrderingExternal::LamellarOrderingExternal()
    : perpDirection_(0),
      interfaceWidth_(),
      externalParameter_(),
      periodicity_(),
      boundaryPtr_(0),
      nAtomType_(0),
      isInitialized_(false)
   {  setClassName("LamellarOrderingExternal"); }

   /*
   * Copy constructor.
   */
   LamellarOrderingExternal::LamellarOrderingExternal(const LamellarOrderingExternal& other)
    : perpDirection_(other.perpDirection_),
      interfaceWidth_(other.interfaceWidth_),
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
   LamellarOrderingExternal& LamellarOrderingExternal::operator = (const LamellarOrderingExternal& other)
   {
      perpDirection_ = other.perpDirection_;
      interfaceWidth_  = other.interfaceWidth_;
      externalParameter_ = other.externalParameter_;
      periodicity_ = other.periodicity_;
      boundaryPtr_ = other.boundaryPtr_;
      nAtomType_ = other.nAtomType_;
      isInitialized_ = other.isInitialized_;
      for (int i=0; i < nAtomType_; ++i) {
        prefactor_[i] = other.prefactor_[i];
      }
      return *this;
   }

   /*
   * Set nAtomType
   */
   void LamellarOrderingExternal::setNAtomType(int nAtomType)
   {
      if (nAtomType <= 0) {
         UTIL_THROW("nAtomType <= 0");
      }
      if (nAtomType > MaxAtomType) {
         UTIL_THROW("nAtomType > LamellarOrderingExternal::MaxAtomType");
      }
      nAtomType_ = nAtomType;
   }

   void LamellarOrderingExternal::setExternalParameter(double externalParameter)
   {
      UTIL_CHECK(isInitialized_);
      externalParameter_ = externalParameter;
   }


   /*
   * Set pointer to the Boundary.
   */
   void LamellarOrderingExternal::setBoundary(Boundary &boundary)
   {  boundaryPtr_ = &boundary; }

   /*
   * Read potential parameters from file.
   */
   void LamellarOrderingExternal::readParameters(std::istream &in)
   {
      if (nAtomType_ <= 0) {
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
      read<double>(in, "interfaceWidth", interfaceWidth_);
      read<int>(in, "periodicity", periodicity_);
      isInitialized_ = true;
   }

   /*
   * Load internal state from an archive.
   */
   void LamellarOrderingExternal::loadParameters(Serializable::IArchive &ar)
   {
      UTIL_CHECK(nAtomType_ > 0);
      loadParameter<int>(ar, "perpDirection", perpDirection_);
      if (perpDirection_ < 0 || perpDirection_ >= Dimension) {
         UTIL_THROW("Invalid index for perpendicular direction.");
      }
      prefactor_.allocate(nAtomType_);
      loadDArray<double>(ar, "prefactor", prefactor_, nAtomType_);
      loadParameter<double>(ar, "externalParameter", externalParameter_);
      loadParameter<double>(ar, "interfaceWidth", interfaceWidth_);
      loadParameter<int>(ar, "periodicity", periodicity_);
      isInitialized_ = true;
   }

   /*
   * Save internal state to an archive.
   */
   void LamellarOrderingExternal::save(Serializable::OArchive &ar)
   {
      UTIL_CHECK(nAtomType_ > 0);
      UTIL_CHECK(isInitialized_);
      ar << perpDirection_;
      ar << prefactor_;
      ar << externalParameter_;
      ar << interfaceWidth_;
      ar << periodicity_;
   }

   /*
   * Set a potential energy parameter, identified by a string.
   */
   void LamellarOrderingExternal::set(std::string name, double value)
   {
      if (name == "externalParameter") {
         externalParameter_ = value;
      } else
      if (name == "interfaceWidth") {
         interfaceWidth_ = value;
      } else
      if (name == "periodicity") {
         periodicity_ = value;
      } else {
         UTIL_THROW("Unrecognized parameter name");
      }
   }

   /*
   * Get a parameter value, identified by a string.
   */
   double LamellarOrderingExternal::get(std::string name) const
   {
      double value = 0.0;
      if (name == "externalParameter") {
         value = externalParameter_;
      } else
      if (name == "interfaceWidth") {
         value = interfaceWidth_;
      } else
      if (name == "periodicity") {
         value = periodicity_;
      } else {
         UTIL_THROW("Unrecognized parameter name");
      }
      return value;
   }

   /*
   * Return external parameter.
   */
   double LamellarOrderingExternal::externalParameter() const
   {  return externalParameter_; }

   /*
   * Return name string "LamellarOrderingExternal".
   */
   std::string LamellarOrderingExternal::className() const
   {  return std::string("LamellarOrderingExternal"); }

}
