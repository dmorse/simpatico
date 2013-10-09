#ifndef INTER_COSINE_SQ_ANGLE_CPP
#define INTER_COSINE_SQ_ANGLE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, Jian Qin and David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "CosineSqAngle.h"
#include <util/math/Constants.h>
#include <util/random/Random.h>

namespace Inter
{

   using namespace Util;

   /* 
   * Constructor.
   */
   CosineSqAngle::CosineSqAngle()
    : nAngleType_(0)
   {
      setClassName("CosineSqAngle");
      for (int i = 0; i < MaxNAngleType; ++i) {
         kappa_[i]     = 0.0;
         theta0_[i]    = 0.0;
         cosTheta0_[i] = 1.0;
      }
   }

   /* 
   * Copy constructor.
   */
   CosineSqAngle::CosineSqAngle(const CosineSqAngle& other)
    : nAngleType_(other.nAngleType_)
   {
      for (int i = 0; i < MaxNAngleType; ++i) {
         kappa_[i]     = other.kappa_[i];
         theta0_[i]    = other.theta0_[i];
         cosTheta0_[i] = other.cosTheta0_[i];
      }
   }

   /* 
   * Assignment.
   */
   CosineSqAngle& CosineSqAngle::operator = (const CosineSqAngle& other)
   {
      nAngleType_ = other.nAngleType_;
      for (int i = 0; i < MaxNAngleType; ++i) {
         kappa_[i]     = other.kappa_[i];
         theta0_[i]    = other.theta0_[i];
         cosTheta0_[i] = other.cosTheta0_[i];
      }
      return *this;
   }

   /* 
   * Set the nAngleType_ member.
   */
   void CosineSqAngle::setNAngleType(int nAngleType)
   {  
      if (nAngleType > MaxNAngleType) {
         UTIL_THROW("nAngleType > CosineSqAngle::MaxNAngleType");
      }
      nAngleType_ = nAngleType;
   }

   /* 
   * Read bend interaction parameters kappa from file.
   */
   void CosineSqAngle::readParameters(std::istream &in) 
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
         cosTheta0_[i] = cos(theta0_[i] / 180.0 * Constants::Pi);
      }
   }

   /*
   * Load internal state from an archive.
   */
   void CosineSqAngle::loadParameters(Serializable::IArchive &ar)
   {
      ar >> nAngleType_; 
      if (nAngleType_ == 0) {
         UTIL_THROW( "nAtomType must be positive");
      }
      // Read parameters
      loadCArray<double> (ar, "kappa", kappa_, nAngleType_);
      loadCArray<double>(ar, "theta0", theta0_, nAngleType_);
      ar.unpack(cosTheta0_, nAngleType_);
   }

   /*
   * Save internal state to an archive.
   */
   void CosineSqAngle::save(Serializable::OArchive &ar)
   {
      ar << nAngleType_;
      ar.pack(kappa_, nAngleType_);
      ar.pack(theta0_, nAngleType_);
      ar.pack(cosTheta0_, nAngleType_);
   }

   /* 
   * Generate a random bond angle chosen from an equilibrium distribution for
   * randomly oriented bonds. 
   *
   * Algorithm: Generate random bond vector and take the angle.
   */
   double CosineSqAngle::randomAngle(
                         Random *random, double beta, int type) const
   {
      double Theta, Theta_min;
      double U, U_min;
      bool chosen = false;
      
      Theta_min = acos(cosTheta0_[type]);
      U_min = energy(cos(Theta_min), type);

      while (chosen == false) {
         Theta = random->uniform(0, Constants::Pi);
         U = energy(cos(Theta), type);
         if (random->uniform() < exp(-beta*(U-U_min))*cos(Theta)/cos(Theta_min))
            chosen = true;
      }

      return Theta;
   }

   /* 
   * Generate a random bond angle cosine chosen from an equilibrium distribution
   * for randomly oriented bonds. 
   *
   * Algorithm: Generate random bond vector and take the angle.
   */
   double CosineSqAngle::randomCosineAngle(
                         Random *random, double beta, int type) const
   {
      double CosineTheta, CosineTheta_min;
      double U, U_min;
      bool chosen = false;
      
      CosineTheta_min = cosTheta0_[type];
      U_min = energy(CosineTheta_min, type);

      while (chosen == false) {
         CosineTheta = random->uniform(-1, 1);
         U = energy(CosineTheta, type);
         if (random->uniform() < exp(-beta*(U-U_min)))
            chosen = true;
      }

      return CosineTheta;
   }

   /*
   * Modify a parameter, identified by a string.
   */
   void CosineSqAngle::set(std::string name, int type, double value)
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
   double CosineSqAngle::get(std::string name, int type) const
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
   * Return name string "CosineSqAngle" for this evaluator class.
   */
   std::string CosineSqAngle::className() const
   {  return std::string("CosineSqAngle"); }


} 
#endif
