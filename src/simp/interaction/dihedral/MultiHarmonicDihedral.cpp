/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "MultiHarmonicDihedral.h"
#include <util/format/Int.h>
#include <util/format/Dbl.h>

#include <iomanip>
namespace Simp
{

   using namespace Util;

   /* 
   * Constructor.
   */
   MultiHarmonicDihedral::MultiHarmonicDihedral()
    : nDihedralType_(0)
   { 
      setClassName("MultiHarmonicDihedral");
   }

   /* 
   * Copy constructor.
   */
   MultiHarmonicDihedral::MultiHarmonicDihedral(const MultiHarmonicDihedral& other)
    : nDihedralType_(other.nDihedralType_)
   {  
      Parameter* ptr;
      const Parameter* otherPtr;
      if (nDihedralType_ > 0) {
         parameters_.allocate(nDihedralType_);
         for (int i = 0; i < nDihedralType_; ++i) {
            ptr = &(parameters_[i]);
            otherPtr = &(other.parameters_[i]);
            ptr->k0 = otherPtr->k0;
            ptr->k1 = otherPtr->k1;
            ptr->k2 = otherPtr->k2;
            ptr->k3 = otherPtr->k3;
            ptr->k4 = otherPtr->k4;
            ptr->init();
         }
      }
   }

   /* 
   * Assignment.
   */
   MultiHarmonicDihedral& 
   MultiHarmonicDihedral::operator = (const MultiHarmonicDihedral& other)
   {
      if (!parameters_.isAllocated()) {
         nDihedralType_ = other.nDihedralType_;
         if (nDihedralType_ > 0) {
            parameters_.allocate(nDihedralType_);
         }
      } else {
         if (nDihedralType_ != other.nDihedralType_) {
           UTIL_THROW("Unequal values of nDihedralType_");
         }
      }
      if (nDihedralType_ > 0) {
         Parameter* ptr;
         const Parameter* otherPtr;
         for (int i = 0; i < nDihedralType_; ++i) {
            ptr = &(parameters_[i]);
            otherPtr = &(other.parameters_[i]);
            ptr->k0 = otherPtr->k0;
            ptr->k1 = otherPtr->k1;
            ptr->k2 = otherPtr->k2;
            ptr->k3 = otherPtr->k3;
            ptr->k4 = otherPtr->k4;
            ptr->init();
         }
      }
      return *this;
   }

   /* 
   * Set the nDihedralType_ member.
   */
   void MultiHarmonicDihedral::setNDihedralType(int nDihedralType)
   {  
      if (nDihedralType <= 0) {
         UTIL_THROW("nDihedralType must be positive");
      }
      nDihedralType_ = nDihedralType;
      parameters_.allocate(nDihedralType_);
   }

   /* 
   * Read Fourier coefficients from file.
   */
   void MultiHarmonicDihedral::readParameters(std::istream &in) 
   {
      if (nDihedralType_ <= 0) {
         UTIL_THROW("nDihedralType must be set before readParam");
      }
      Parameter* ptr;
      int j;
      for (int i = 0; i < nDihedralType_; ++i) {
         in >> j;
         if (i != j) UTIL_THROW("Inconsistent dihedral type index");
         ptr = &(parameters_[i]);
         in >> ptr->k0 >> ptr->k1 >> ptr->k2 
            >> ptr->k3 >> ptr->k4;
         ptr->init();
      }
   }

   /* 
   * Write Fourier coefficients to file.
   */
   void MultiHarmonicDihedral::writeParam(std::ostream &out) 
   {
      Parameter* ptr;
      out.setf(std::ios_base::fixed);
      for (int i = 0; i < nDihedralType_; ++i) {
         ptr = &parameters_[i];
         out << Int(i, 5) << "  " << std::setprecision(6)
             << std::setw(10) << ptr->k0   << "  "
             << std::setw(10) << ptr->k1   << "  "
             << std::setw(10) << ptr->k2   << "  "
             << std::setw(10) << ptr->k3   << "  "
             << std::setw(10) << ptr->k4   << "  "
             << std::endl;
      }
   }

   /*
   * Load internal state from an archive.
   */
   void MultiHarmonicDihedral::loadParameters(Serializable::IArchive &ar)
   {
      ar >> nDihedralType_; 
      if (nDihedralType_ <= 0) {
         UTIL_THROW( "nDihedralType must be positive");
      }
      parameters_.allocate(nDihedralType_);
      Parameter* ptr;
      for (int i = 0; i < nDihedralType_; ++i) {
         ptr = &(parameters_[i]);
         ar >> ptr->k0;
         ar >> ptr->k1;
         ar >> ptr->k2;
         ar >> ptr->k3;
         ar >> ptr->k4;
         ar >> ptr->a0;
         ar >> ptr->a1;
         ar >> ptr->a2;
         ar >> ptr->a3;
         ar >> ptr->a4;
         ar >> ptr->g2;
         ar >> ptr->g3;
         ar >> ptr->g4;
      }
   }

   /*
   * Save internal state to an archive.
   */
   void MultiHarmonicDihedral::save(Serializable::OArchive &ar)
   {
      ar << nDihedralType_;
      Parameter* ptr;
      for (int i = 0; i < nDihedralType_; ++i) {
         ptr = &(parameters_[i]);
         ar << ptr->k0;
         ar << ptr->k1;
         ar << ptr->k2;
         ar << ptr->k3;
         ar << ptr->k4;
         ar << ptr->a0;
         ar << ptr->a1;
         ar << ptr->a2;
         ar << ptr->a3;
         ar << ptr->a4;
         ar << ptr->g2;
         ar << ptr->g3;
         ar << ptr->g4;
      }
   }

   /*
   * Modify a parameter, identified by a string.
   */
   void MultiHarmonicDihedral::set(std::string name, int type, double value)
   {
      Parameter* ptr = &(parameters_[type]);
      if (name == "k0") {
        ptr->k0 = value;
      } else 
      if (name == "k1") {
        ptr->k1 = value;
      } else 
      if (name == "k2") {
        ptr->k2 = value;
      } else 
      if (name == "k3") {
        ptr->k3 = value;
      } else 
      if (name == "k4") {
        ptr->k4 = value;
      } else {
         UTIL_THROW("Unrecognized parameter name");
      }
      ptr->init();
   }

   /*
   * Get a parameter value, identified by a string.
   */
   double MultiHarmonicDihedral::get(std::string name, int type) const
   {
      double value = 0.0;
      const Parameter* ptr = &(parameters_[type]);
      if (name == "k0") {
         value = ptr->k0;
      } else 
      if (name == "k1") {
         value = ptr->k1;
      } else 
      if (name == "k2") {
         value = ptr->k2;
      } else 
      if (name == "k3") {
         value = ptr->k3;
      } else 
      if (name == "k4") {
         value = ptr->k4;
      } else {
         UTIL_THROW("Unrecognized parameter name");
      }
      return value;
   }

   /*
   * Return class name string "MultiHarmonicDihedral".
   */
   std::string MultiHarmonicDihedral::className() const
   {  return std::string("MultiHarmonicDihedral"); }

} 
