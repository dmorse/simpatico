#ifndef HARMONIC_L0_BOND_CPP
#define HARMONIC_L0_BOND_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
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
   void HarmonicL0Bond::readParam(std::istream &in) 
   {

      // Precondition
      if (nBondType_ <= 0) {
         UTIL_THROW("nBondType must be set before readParam");
      }

      // Read parameters
      //readBegin(in, "HarmonicL0Bond");
      readCArray<double>(in, "kappa",  kappa_,  nBondType_);
      //readEnd(in);
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
   * Return name string "HarmonicL0Bond" for this evaluator class.
   */
   std::string HarmonicL0Bond::className() const
   {  return std::string("HarmonicL0Bond"); }


} 
#endif
