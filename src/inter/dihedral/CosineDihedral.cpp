#ifndef INTER_COSINE_DIHEDRAL_CPP
#define INTER_COSINE_DIHEDRAL_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, Jian Qin and David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "CosineDihedral.h"

namespace Inter
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
      // Preconditions
      if (nDihedralType_ <= 0) {
         UTIL_THROW("nDihedralType must be set before readParam");
      }

      // Read parameters
      //readBegin(in, "CosineDihedral");
      readCArray<double>(in, "kappa",  kappa_,  nDihedralType_);
      //readEnd(in);
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
#endif
