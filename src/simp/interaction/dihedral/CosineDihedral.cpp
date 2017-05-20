/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "CosineDihedral.h"

namespace Simp
{

   using namespace Util;

   /* 
   * Constructor.
   */
   CosineDihedral::CosineDihedral()
    : nDihedralType_(0)
   { 
      setClassName("CosineDihedral");
      for (int i = 0; i < MaxNDihedralType; ++i) {
         kappa_[i] = 0.0; 
      }
   }

   /* 
   * Copy constructor.
   */
   CosineDihedral::CosineDihedral(const CosineDihedral& other)
    : nDihedralType_(other.nDihedralType_)
   { for (int i = 0; i < nDihedralType_; ++i) kappa_[i] = other.kappa_[i]; }

   /* 
   * Assignment.
   */
   CosineDihedral& CosineDihedral::operator = (const CosineDihedral& other)
   {
      nDihedralType_ = other.nDihedralType_;
      for (int i = 0; i < nDihedralType_; ++i) kappa_[i] = other.kappa_[i];
      return *this;
   }

   /* 
   * Set the nDihedralType_ member.
   */
   void CosineDihedral::setNDihedralType(int nDihedralType)
   {  
      if (nDihedralType > MaxNDihedralType) {
         UTIL_THROW("nDihedralType > CosineDihedral::MaxNDihedralType");
      }
      nDihedralType_ = nDihedralType;
   }

   /* 
   * Read bend interaction parameters kappa from file.
   */
   void CosineDihedral::readParameters(std::istream &in) 
   {
      if (nDihedralType_ <= 0) {
         UTIL_THROW("nDihedralType must be set before readParam");
      }
      readCArray<double>(in, "kappa",  kappa_,  nDihedralType_);
   }

   /*
   * Load internal state from an archive.
   */
   void CosineDihedral::loadParameters(Serializable::IArchive &ar)
   {
      ar >> nDihedralType_; 
      if (nDihedralType_ <= 0) {
         UTIL_THROW( "nDihedralType must be positive");
      }
      loadCArray<double> (ar, "kappa", kappa_, nDihedralType_);
   }

   /*
   * Save internal state to an archive.
   */
   void CosineDihedral::save(Serializable::OArchive &ar)
   {
      ar << nDihedralType_;
      ar.pack(kappa_, nDihedralType_);
   }

   /*
   * Modify a parameter, identified by a string.
   */
   void CosineDihedral::set(std::string name, int type, double value)
   {
      if (name == "kappa") {
         kappa_[type] = value;
      } else {
         UTIL_THROW("Unrecognized parameter name");
      }
   }

   /*
   * Get a parameter value, identified by a string.
   */
   double CosineDihedral::get(std::string name, int type) const
   {
      double value = 0.0;
      if (name == "kappa") {
         value = kappa_[type];
      } else {
         UTIL_THROW("Unrecognized parameter name");
      }
      return value;
   }

   /*
   * Return name string "CosineDihedral" for this evaluator class.
   */
   std::string CosineDihedral::className() const
   {  return std::string("CosineDihedral"); }

} 
