/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "HarmonicBond.h"
#include <util/random/Random.h>

namespace Simp
{

   using namespace Util;

   /* 
   * Constructor.
   */
   HarmonicBond::HarmonicBond()
    : nBondType_(0)
   {
      setClassName("HarmonicBond"); 
      for (int i = 0; i < MaxNBondType; ++i) {
         kappa_[i]  =  0.0;
         length_[i] =  0.0;
      }
   }

   /* 
   * Copy constructor.
   */
   HarmonicBond::HarmonicBond(const HarmonicBond& other)
    : nBondType_(other.nBondType_)
   {
      for (int i = 0; i < nBondType_; ++i) {
         kappa_[i]  =  other.kappa_[i];
         length_[i] =  other.length_[i];
      }
   }

   /* 
   * Assignment.
   */
   HarmonicBond& HarmonicBond::operator = (const HarmonicBond& other)
   {
      nBondType_  = other.nBondType_;
      for (int i = 0; i < nBondType_; ++i) {
         kappa_[i]  = other.kappa_[i];
         length_[i] = other.length_[i];
      }
      return *this;
   }

   /* 
   * Set the nBondType_ member
   */
   void HarmonicBond::setNBondType(int nBondType)
   {  
      if (nBondType > MaxNBondType) {
         UTIL_THROW("nBondType > HarmonicBond::MaxNBondType");
      }
      nBondType_ = nBondType;
   }

   /* 
   * Read bond interaction parameters kappa and length from file
   */
   void HarmonicBond::readParameters(std::istream &in) 
   {
      UTIL_CHECK(nBondType_ > 0);
      readCArray<double>(in, "kappa",  kappa_,  nBondType_);
      readCArray<double>(in, "length", length_, nBondType_);
   }
   
   /*
   * Load internal state from an archive.
   */
   void HarmonicBond::loadParameters(Serializable::IArchive &ar)
   {
      UTIL_CHECK(nBondType_ > 0);
      loadCArray<double>(ar, "kappa", kappa_, nBondType_);
      loadCArray<double>(ar, "length", length_, nBondType_);
   }

   /*
   * Save internal state to an archive.
   */
   void HarmonicBond::save(Serializable::OArchive &ar)
   {
      UTIL_CHECK(nBondType_ > 0);
      ar.pack(kappa_, nBondType_);
      ar.pack(length_, nBondType_);
   }

   /* 
   * Generate a random bond length chosen from an equilibrium distribution for
   * randomly oriented bonds. 
   *
   * The desired probability distribution is proportional to x*x*G(x) for x > 0,
   * where G(x) is a Gaussian distribution with a mean x0=length_[type] and a 
   * standard deviation sd = sqrt(temperature/kappa_[type]).
   *
   * Algorithm:
   *   - Generate trial bond length x from the Gaussian distribution G(x).
   *   - If x < 0, generate another trial bond length.
   *   - Let y = x*x/xm*xm, where xm = x0 + 7.0*sd is a maximum length.
   *     The algorithm is correct only for y < 1, or x < xm. The distribution
   *     G(y) for generated values of y is proportional to G(x)/x.
   *   - Accept y if it is greater than a uniform random number in [0,1]. 
   *     For y < 1, or x < xm, the probability of accepting y is equal to y. 
   *     The distribution A(y) of accepted values of y for y < 1 is thus 
   *     proportional to y*G(y), or to x*G(x).
   *   - If y is accepted, return x. The distribution A(x) for accepted 
   *     values of x is proportional to x*x*G(x) for x < xm.
   *
   * Limitations:
   *     The distribution is correct only for x < xm. The probability of
   *     generating a value x >= xm is, however, less than exp(-49/2).
   */
   double HarmonicBond::randomBondLength(
                         Random *random, double beta, int type) const
   {
   
      double sd  = 1.0/sqrt( beta*kappa_[type] );
      double x0  = length_[type];
      double xm2 = length_[type] + 7.0*sd;
             xm2 = xm2*xm2;
      double x;
   
      // Generate trials until one is accepted
      for (;;) { 
         x = x0 + sd*random->gaussian();
         if (x < 0) continue;
         if ( x*x/xm2 >= random->uniform(0.0,1.0) ) {
            return x;
         };
      }
   }

   /*
   * Modify a parameter, identified by a string.
   */
   void HarmonicBond::set(std::string name, int type, double value)
   {
      if (name == "kappa") {
         kappa_[type] = value;
      } else
      if (name == "length") {
         length_[type] = value;
      } else {
         UTIL_THROW("Unrecognized parameter name");
      }
   }

   /*
   * Get a parameter value, identified by a string.
   */
   double HarmonicBond::get(std::string name, int type) const
   {
      double value = 0.0;
      if (name == "kappa") {
         value = kappa_[type];
      } else
      if (name == "length") {
         value = length_[type];
      } else {
         UTIL_THROW("Unrecognized parameter name");
      }
      return value;
   }

} 
