/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "HarmonicAngle.h"
#include <util/math/Constants.h>
#include <util/random/Random.h>


namespace Simp
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
      // Precondition
      UTIL_CHECK(nAngleType_ > 0);

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
      UTIL_CHECK(nAngleType_ > 0);
      loadCArray<double> (ar, "kappa", kappa_, nAngleType_);
      loadCArray<double>(ar, "theta0", theta0_, nAngleType_);
   }

   /*
   * Save internal state to an archive.
   */
   void HarmonicAngle::save(Serializable::OArchive &ar)
   {
      UTIL_CHECK(nAngleType_ > 0);
      ar.pack(kappa_, nAngleType_);
      ar.pack(theta0_, nAngleType_);
   }

   /* 
   * Generate a random bond angle chosen from an equilibrium distribution for
   * randomly oriented bonds. 
   *
   * Algorithm: Generate random bond vector and take the angle.
   */
   double HarmonicAngle::randomAngle(
                         Random *random, double beta, int type) const
   {
      double Theta, Theta_min;
      double U, U_min;
      bool chosen = false;
      
      Theta_min = theta0_[type];
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
   double HarmonicAngle::randomCosineAngle(
                         Random *random, double beta, int type) const
   {
      double CosineTheta, CosineTheta_min;
      double U, U_min;
      bool chosen = false;
      
      CosineTheta_min = cos(theta0_[type]);
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
