#ifdef  MCMD_ANGLE
#ifndef COSINE_ANGLE_CPP
#define COSINE_ANGLE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, Jian Qin and David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "CosineAngle.h"

namespace McMd
{

   using namespace Util;

   /* 
   * Constructor.
   */
   CosineAngle::CosineAngle()
    : nAngleType_(0)
   { 
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
   void CosineAngle::readParam(std::istream &in) 
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
   * Return name string "CosineAngle" for this evaluator class.
   */
   std::string CosineAngle::className() const
   {  return std::string("CosineAngle"); }

}
 
#endif
#endif
