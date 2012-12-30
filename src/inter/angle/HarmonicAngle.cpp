#ifndef HARMONIC_ANGLE_CPP
#define HARMONIC_ANGLE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, Jian Qin and David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "HarmonicAngle.h"
#include <util/math/Constants.h>

namespace Inter
{

   using namespace Util;

   /* 
   * Constructor.
   */
   HarmonicAngle::HarmonicAngle()
    : nAngleType_(0)
   {
      setClassName("HarmonicAngle");
      for (int i = 0; i < MaxNAngleType; ++i) {
         kappa_[i]     = 0.0;
         theta0_[i]    = 0.0;
      }
   }

   /* 
   * Copy constructor.
   */
   HarmonicAngle::HarmonicAngle(const HarmonicAngle& other)
    : nAngleType_(other.nAngleType_)
   {
      for (int i = 0; i < MaxNAngleType; ++i) {
         kappa_[i]     = other.kappa_[i];
         theta0_[i]    = other.theta0_[i];
      }
   }

   /* 
   * Assignment.
   */
   HarmonicAngle& HarmonicAngle::operator = (const HarmonicAngle& other)
   {
      nAngleType_ = other.nAngleType_;
      for (int i = 0; i < MaxNAngleType; ++i) {
         kappa_[i]     = other.kappa_[i];
         theta0_[i]    = other.theta0_[i];
      }
      return *this;
   }

   /* 
   * Set the nAngleType_ member.
   */
   void HarmonicAngle::setNAngleType(int nAngleType)
   {  
      if (nAngleType > MaxNAngleType) {
         UTIL_THROW("nAngleType > HarmonicAngle::MaxNAngleType");
      }
      nAngleType_ = nAngleType;
   }

   /* 
   * Read bend interaction parameters kappa from file.
   */
   void HarmonicAngle::readParameters(std::istream &in) 
   {
      // Preconditions
      if (nAngleType_ <= 0) {
         UTIL_THROW("nAngleType must be set before readParam");
      }

      // Read parameters
      readCArray<double>(in, "kappa",  kappa_,  nAngleType_);
      readCArray<double>(in, "theta0", theta0_, nAngleType_);

      // Convert from degrees to radians.
      for (int i = 0; i < nAngleType_; ++i) {
         theta0_[i] = theta0_[i] * Constants::Pi / 180.0;
      }
   }

   /*
   * Load internal state from an archive.
   */
   void HarmonicAngle::loadParameters(Serializable::IArchive &ar)
   {
      ar >> nAngleType_; 
      if (nAngleType_ == 0) {
         UTIL_THROW( "nAtomType must be positive");
      }
      // Read parameters
      loadCArray<double> (ar, "kappa", kappa_, nAngleType_);
      loadCArray<double>(ar, "theta0", theta0_, nAngleType_);
   }

   /*
   * Save internal state to an archive.
   */
   void HarmonicAngle::save(Serializable::OArchive &ar)
   {
      ar << nAngleType_;
      ar.pack(kappa_, nAngleType_);
      ar.pack(theta0_, nAngleType_);
   }

   /*
   * Modify a parameter, identified by a string.
   */
   void HarmonicAngle::set(std::string name, int type, double value)
   {
      if (name == "kappa") {
         kappa_[type] = value;
      } else
      if (name == "theta0") {
         theta0_[type] = value;
      } else {
         UTIL_THROW("Unrecognized parameter name");
      }
   }

   /*
   * Get a parameter value, identified by a string.
   */
   double HarmonicAngle::get(std::string name, int type) const
   {
      double value = 0.0;
      if (name == "kappa") {
         value = kappa_[type];
      } else
      if (name == "theta0") {
         value = theta0_[type];
      } else {
         UTIL_THROW("Unrecognized parameter name");
      }
      return value;
   }

   /*
   * Return name string "HarmonicAngle" for this evaluator class.
   */
   std::string HarmonicAngle::className() const
   {  return std::string("HarmonicAngle"); }


} 
#endif
