#ifndef INTER_FENE_BOND_CPP
#define INTER_FENE_BOND_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "FeneBond.h"
#include <util/random/Random.h>

namespace Inter
{

   using namespace Util;

   /*
   * Constructor.
   */
   FeneBond::FeneBond()
    : nBondType_(0)
   {
      setClassName("FeneBond"); 
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
   * Read bond interaction parameters kappa and r0 from file
   */
   void FeneBond::readParameters(std::istream &in)
   {

      // Preconditions
      if (nBondType_ <= 0) {
         UTIL_THROW("nBondType must be set before readParam");
      }

      // Read parameters
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
         //std::cout << "Distance cutoff = " << rl_[i] << std::endl;
         //std::cout << "Energy cutoff   = " << energyCutoff_[i] << std::endl;
      }
   }

   /*
   * Load internal state from an archive.
   */
   void FeneBond::loadParameters(Serializable::IArchive &ar)
   {
      ar >> nBondType_; 
      if (nBondType_ == 0) {
         UTIL_THROW( "nAtomType must be positive");
      }
      // Read parameters
      loadCArray<double> (ar, "kappa", kappa_, nBondType_);
      loadCArray<double>(ar, "r0", r0_, nBondType_);
      ar >> forceCutoff_;
      ar.unpack(r0Sq_, nBondType_);
      ar.unpack(r0SqInv_, nBondType_);
      ar.unpack(ce_, nBondType_);
      ar.unpack(rl_, nBondType_);
      ar.unpack(rlSq_, nBondType_);
      ar.unpack(energyCutoff_, nBondType_);
   }

   /*
   * Save internal state to an archive.
   */
   void FeneBond::save(Serializable::OArchive &ar)
   {
      ar << nBondType_;
      ar.pack(kappa_, nBondType_);
      ar.pack(r0_, nBondType_);
      ar << forceCutoff_;
      ar.pack(r0Sq_, nBondType_);
      ar.pack(r0SqInv_, nBondType_);
      ar.pack(ce_, nBondType_);
      ar.pack(rl_, nBondType_);
      ar.pack(rlSq_, nBondType_);
      ar.pack(energyCutoff_, nBondType_);
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
         if (random->uniform(0.0,1.0) < ratio) {
            return sqrt(rSq);
         };
      }
   }

   /*
   * Modify a parameter, identified by a string.
   */
   void FeneBond::set(std::string name, int type, double value)
   {
      if (name == "kappa") {
         kappa_[type] = value;
      } else
      if (name == "r0") {
         r0_[type] = value;
      } else {
         UTIL_THROW("Unrecognized parameter name");
      }
      r0Sq_[type] = r0_[type]*r0_[type];
      r0SqInv_[type] = 1.0/r0Sq_[type];
      ce_[type] = -0.5*r0Sq_[type]*kappa_[type];
      //a_[type] = exp(-2.0*energyCutoff_/(r0Sq_[type]*kappa_[type]));
      //rl_[type] = r0_[type]*sqrt(1-a_[type]);
      double y = 0.5*kappa_[type]*r0_[type]/forceCutoff_;
      rl_[type] = r0_[type]*(sqrt(1.0 + y*y) - y);
      rlSq_[type] = rl_[type]*rl_[type];
      double g = 1.0 - rlSq_[type]*r0SqInv_[type];
      energyCutoff_[type] = ce_[type]*log(g);
      //std::cout << "Force cutoff    = " << kappa_[type]*rl_[type]/g << std::endl;
      //std::cout << "Distance cutoff = " << rl_[type] << std::endl;
      //std::cout << "Energy cutoff   = " << energyCutoff_[type] << std::endl;
   }

   /*
   * Get a parameter value, identified by a string.
   */
   double FeneBond::get(std::string name, int type) const
   {
      double value;
      if (name == "kappa") {
         value = kappa_[type];
      } else
      if (name == "r0") {
         value = r0_[type];
      } else {
         UTIL_THROW("Unrecognized parameter name");
      }
      return value;
   }

}
#endif

