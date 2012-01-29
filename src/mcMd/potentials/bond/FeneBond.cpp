#ifndef FENE_BOND_CPP
#define FENE_BOND_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "FeneBond.h"
#include <util/random/Random.h>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   FeneBond::FeneBond()
    : nBondType_(0)
   {
      for (int i = 0; i < MaxNBondType; ++i) {
         kappa_[i]   =  0.0;
         r0_[i]      =  0.0;
         r0Sq_[i]    =  0.0;
         r0SqInv_[i] =  0.0;
         ce_[i]      =  0.0;
         rlSq_[i]    =  0.0;
         rl_[i]      =  0.0;
         energyCutoff_[i] =  0.0;
         //a_[i]       =  0.0;
      }
   }

   /*
   * Copy constructor.
   */
   FeneBond::FeneBond(const FeneBond& other)
    : nBondType_(other.nBondType_)
   {
      for (int i = 0; i < nBondType_; ++i) {
         kappa_[i]   =  other.kappa_[i];
         r0_[i]      =  other.r0_[i];
         r0Sq_[i]    =  other.r0Sq_[i];
         r0SqInv_[i] =  other.r0SqInv_[i];
         ce_[i]      =  other.ce_[i];
         rlSq_[i]    =  other.rlSq_[i];
         //a_[i]     =  other.a_[i];
         energyCutoff_[i] = other.energyCutoff_[i];
         rl_[i]      =  other.rl_[i];
      }
   }

   /*
   * Assignment.
   */
   FeneBond& FeneBond::operator = (const FeneBond& other)
   {
      nBondType_   = other.nBondType_;
      forceCutoff_ = other.forceCutoff_;
      for (int i = 0; i < nBondType_; ++i) {
         kappa_[i]   = other.kappa_[i];
         r0_[i]      = other.r0_[i];
         r0Sq_[i]    = other.r0Sq_[i];
         r0SqInv_[i] =  other.r0SqInv_[i];
         ce_[i]      =  other.ce_[i];
         rlSq_[i]    =  other.rlSq_[i];    
         rl_[i]      =  other.rl_[i];
         energyCutoff_[i] = other.energyCutoff_[i];
         //a_[i]       =  other.a_[i];
      }
      return *this;
   }

   /*
   * Set the nBondType_ member
   */
   void FeneBond::setNBondType(int nBondType)
   {
      if (nBondType > MaxNBondType) {
         UTIL_THROW("nBondType > FeneBond::MaxNBondType");
      }
      nBondType_ = nBondType;
   }

   /*
   * Read bond interaction parameters kappa and length from file
   */
   void FeneBond::readParam(std::istream &in)
   {

      // Preconditions
      if (nBondType_ <= 0) {
         UTIL_THROW("nBondType must be set before readParam");
      }

      // Read parameters
      //readBegin(in, "FeneBond");
      readCArray<double>(in, "kappa",  kappa_,  nBondType_);
      readCArray<double>(in, "r0", r0_, nBondType_);
      //read<double>(in, "energyCutoff", energyCutoff_);
      read<double>(in, "forceCutoff", forceCutoff_);

      double y, g;
      for (int i = 0; i < nBondType_; ++i) {
         r0Sq_[i] = r0_[i]*r0_[i];
         r0SqInv_[i] = 1.0/r0Sq_[i];
         ce_[i] = -0.5*r0Sq_[i]*kappa_[i];
         //a_[i] = exp(-2.0*energyCutoff_/(r0Sq_[i]*kappa_[i]));
         //rl_[i] = r0_[i]*sqrt(1-a_[i]);
         y = 0.5*kappa_[i]*r0_[i]/forceCutoff_;
         rl_[i] = r0_[i]*(sqrt(1.0 + y*y) - y);
         rlSq_[i] = rl_[i]*rl_[i];
         g = 1.0 - rlSq_[i]*r0SqInv_[i];
         energyCutoff_[i] = ce_[i]*log(g);
         //std::cout << "Force cutoff    = " << kappa_[i]*rl_[i]/g << std::endl;
         std::cout << "Distance cutoff = " << rl_[i] << std::endl;
         std::cout << "Energy cutoff   = " << energyCutoff_[i] << std::endl;
      }

      //readEnd(in);
   }

   /*
   * Generate a random bond length chosen from an equilibrium distribution for
   * randomly oriented bonds.
   *
   * Algorithm:
   *
   *   While trial not accepted:
   *
   *      -Choose random trial vector of squared length rSq from the equilibrium
   *       distribution for a harmonic spring with the same spring constant.
   *      -Reject if rSq > r0Sq
   *      -Accept with probability exp[-beta[U_Fene(r) - U_Harmonic(r)].
   */
   double FeneBond::randomBondLength(
                         Random *random, double beta, int type) const
   {

      double x, y, z, sigSq, rSq, eHarm, eFene, ratio;
      sigSq = 1.0/(beta*kappa_[type]);

      // Generate trials until one is accepted
      for (;;) {
         x = random->gaussian();
         y = random->gaussian();
         z = random->gaussian();
         rSq = sqrt(x*x + y*y + z*z)*sigSq;
         if (rSq >= r0Sq_[type]) continue;
         eHarm = 0.5*kappa_[type]*rSq;
         eFene = energy(rSq, type);
         ratio = exp(-beta*(eFene - eHarm));
         if (random->getFloat(0.0,1.0) < ratio) {
            return sqrt(rSq);
         };
      }
   }

   /*
   * Return name string "FeneBond" for this evaluator class.
   */
   std::string FeneBond::className() const
   {  return std::string("FeneBond"); }

}
#endif

