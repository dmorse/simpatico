#ifdef  MCMD_TETHER
/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "HarmonicTether.h"

namespace Simp
{

   using namespace Util;

   /* 
   * Default constructor.
   */
   HarmonicTether::HarmonicTether()
    : nTetherType_(0)
   {
      for (int i=0; i < MaxNTetherType; ++i) {
         kappa_[i] = 0.0;
      }
   }

   /* 
   * Copy constructor.
   */
   HarmonicTether::HarmonicTether(const HarmonicTether& other)
    : nTetherType_(other.nTetherType_)
   {
      for (int i = 0; i < nTetherType_; ++i) {
         kappa_[i] = other.kappa_[i];
      }
   }

   /* 
   * Assignment.
   */
   HarmonicTether& HarmonicTether::operator = (const HarmonicTether& other)
   {
      nTetherType_ = other.nTetherType_;
      for (int i = 0; i < nTetherType_; ++i) {
         kappa_[i] = other.kappa_[i];
      }
      return *this;
   }

   /* 
   * Read bond interaction parameters from file
   */
   void HarmonicTether::readParameters(std::istream &in) 
   {
      //readBegin(in, "HarmonicTether");
      read<int>(in, "nTetherType", nTetherType_);
      if (nTetherType_ > MaxNTetherType) {
         UTIL_THROW("nTetherType_ > MaxNTetherType");
      }
      readCArray<double>(in, "kappa",  kappa_,  nTetherType_);
      //readEnd(in);
   }

} 
#endif
