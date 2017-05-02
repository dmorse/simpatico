/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>
#include "CosineAngle.h"
#include <util/math/Constants.h>
#include <util/random/Random.h>

namespace Simp
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
      UTIL_CHECK(nAngleType_ > 0);
      readCArray<double>(in, "kappa",  kappa_,  nAngleType_);
   }

   /*
   * Load internal state from an archive.
   */
   void CosineAngle::loadParameters(Serializable::IArchive &ar)
   {
      UTIL_CHECK(nAngleType_ > 0);
      loadCArray<double> (ar, "kappa", kappa_, nAngleType_);
   }

   /*
   * Save internal state to an archive.
   */
   void CosineAngle::save(Serializable::OArchive &ar)
   {  
      UTIL_CHECK(nAngleType_ > 0);
      ar.pack(kappa_, nAngleType_); 
   }

   /* 
   * Generate a random bond angle chosen from an equilibrium distribution for
   * randomly oriented bonds. 
   *
   * Algorithm: Generate random bond vector and take the angle.
   */
   double CosineAngle::randomAngle(Random *random, double beta, int type) const
   {
      double Theta, Theta_min;
      double U, U_min;
      bool chosen = false;
      
      Theta_min = 0;
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
   double CosineAngle::randomCosineAngle(Random *random, double beta, int type) const
   {
      double CosineTheta, CosineTheta_min;
      double U, U_min;
      bool chosen = false;
      
      CosineTheta_min = 1;
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
