#ifndef INTER_HARMONIC_L0_BOND_CPP
#define INTER_HARMONIC_L0_BOND_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "HarmonicL0Bond.h"
#include <util/random/Random.h>

namespace Inter
{

   using namespace Util;

   /* 
   * Constructor.
   */
   HarmonicL0Bond::HarmonicL0Bond()
    : nBondType_(0)
   {
      setClassName("HarmonicL0Bond"); 
      for (int i = 0; i < MaxNBondType; ++i) {
         kappa_[i]  =  0.0;
      }
   }

   /* 
   * Copy constructor.
   */
   HarmonicL0Bond::HarmonicL0Bond(const HarmonicL0Bond& other)
    : nBondType_(other.nBondType_)
   {
      for (int i = 0; i < nBondType_; ++i) {
         kappa_[i]  =  other.kappa_[i];
      }
   }

   /* 
   * Assignment.
   */
   HarmonicL0Bond& HarmonicL0Bond::operator = (const HarmonicL0Bond& other)
   {
      nBondType_   = other.nBondType_;
      for (int i = 0; i < nBondType_; ++i) {
         kappa_[i]  = other.kappa_[i];
      }
      return *this;
   }

   /* 
   * Set the nBondType_ member
   */
   void HarmonicL0Bond::setNBondType(int nBondType)
   {  
      if (nBondType > MaxNBondType) {
         UTIL_THROW("nBondType > HarmonicL0Bond::MaxNBondType");
      }
      nBondType_ = nBondType;
   }

   /* 
   * Read bond interaction parameters kappa and length from file
   */
   void HarmonicL0Bond::readParameters(std::istream &in) 
   {

      // Precondition
      if (nBondType_ <= 0) {
         UTIL_THROW("nBondType must be set before readParam");
      }
      readCArray<double>(in, "kappa",  kappa_,  nBondType_);
   }

   /*
   * Load internal state from an archive.
   */
   void HarmonicL0Bond::loadParameters(Serializable::IArchive &ar)
   {
      ar >> nBondType_; 
      if (nBondType_ == 0) {
         UTIL_THROW( "nBondType must be positive");
      }
      loadCArray<double> (ar, "kappa", kappa_, nBondType_);
   }

   /*
   * Save internal state to an archive.
   */
   void HarmonicL0Bond::save(Serializable::OArchive &ar)
   {
      ar << nBondType_;
      ar.pack(kappa_, nBondType_);
   }

   
   /* 
   * Generate a random bond length chosen from an equilibrium distribution for
   * randomly oriented bonds. 
   *
   * Algorithm: Generate random bond vector and take the sqrt.
   */
   double HarmonicL0Bond::randomBondLength(
                         Random *random, double beta, int type) const
   {
      double x, y, z;
      x = random->gaussian();
      y = random->gaussian();
      z = random->gaussian();
      return sqrt(x*x + y*y + z*z)/sqrt(beta*kappa_[type]);
   }

   /*
   * Modify a parameter, identified by a string.
   */
   void HarmonicL0Bond::set(std::string name, int type, double value)
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
   double HarmonicL0Bond::get(std::string name, int type) const
   {
      double value;
      if (name == "kappa") {
         value = kappa_[type];
      } else {
         UTIL_THROW("Unrecognized parameter name");
      }
      return value;
   }

} 
#endif
