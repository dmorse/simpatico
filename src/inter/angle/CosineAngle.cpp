#ifndef INTER_COSINE_ANGLE_CPP
#define INTER_COSINE_ANGLE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, Jian Qin and David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "CosineAngle.h"

namespace Inter
{

   using namespace Util;

   /* 
   * Constructor.
   */
   CosineAngle::CosineAngle()
    : nAngleType_(0)
   { 
      setClassName("CosineAngle");
      for (int i = 0; i < MaxNAngleType; ++i) {
         kappa_[i] = 0.0; 
      }
   }

   /* 
   * Copy constructor.
   */
   CosineAngle::CosineAngle(const CosineAngle& other)
    : nAngleType_(other.nAngleType_)
   { 
      assert(other.nAngleType_ > 0);
      for (int i = 0; i < nAngleType_; ++i) {
         kappa_[i] = other.kappa_[i]; 
      }
   }

   /* 
   * Destructor.
   */
   CosineAngle::~CosineAngle()
   {}
 
   /* 
   * Assignment.
   */
   CosineAngle& CosineAngle::operator = (const CosineAngle& other)
   {
      assert(other.nAngleType_ > 0);

      nAngleType_ = other.nAngleType_;
      for (int i = 0; i < nAngleType_; ++i) {
         kappa_[i] = other.kappa_[i];
      }
      return *this;
   }

   /* 
   * Set the nAngleType_ member.
   */
   void CosineAngle::setNAngleType(int nAngleType)
   {  
      if (nAngleType > MaxNAngleType) {
         UTIL_THROW("nAngleType > CosineAngle::MaxNAngleType");
      }
      nAngleType_ = nAngleType;
   }

   /* 
   * Read bend interaction parameters kappa from file.
   */
   void CosineAngle::readParameters(std::istream &in) 
   {
      // Preconditions
      if (nAngleType_ <= 0) {
         UTIL_THROW("nAngleType must be set before readParam");
      }

      // Read parameters
      //readBegin(in, "CosineAngle");
      readCArray<double>(in, "kappa",  kappa_,  nAngleType_);
      //readEnd(in);
   }

   /*
   * Modify a parameter, identified by a string.
   */
   void CosineAngle::set(std::string name, int type, double value)
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
   double CosineAngle::get(std::string name, int type) const
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
   * Return name string "CosineAngle" for this evaluator class.
   */
   std::string CosineAngle::className() const
   {  return std::string("CosineAngle"); }

}
 
#endif
