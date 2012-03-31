#ifndef INTER_NOPAIR
#ifndef HOOMD_LJ_SHIFTED_FORCE_PAIR_CPP
#define HOOMD_LJ_SHIFTED_FORCE_PAIR_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "HoomdLJShiftedForcePair.h"


namespace McMd
{
   char classNameHoomdLJShiftedForce[] = "HoomdLJShiftedForcePair";

   /**
   * Default constructor.
   */
   HoomdLJShiftedForcePair::HoomdLJShiftedForcePair()   
    : HoomdPair< EvaluatorPairLJ, gpu_compute_ljtemp_forces,
         classNameHoomdLJShiftedForce >()
   {
   }

   /**
   * Copy constructor
   */
   HoomdLJShiftedForcePair::HoomdLJShiftedForcePair(const HoomdLJShiftedForcePair& other)
    : HoomdPair< EvaluatorPairLJ, gpu_compute_ljtemp_forces,
         classNameHoomdLJShiftedForce >(other)
   {
      for (int i = 0; i < nAtomType_; ++i) {
         for (int j = 0; j < nAtomType_; ++j) {
            epsilon_[i][j] = other.epsilon_[i][j];
            sigma_[i][j] = other.sigma_[i][j];
         }
      }
   }

   /**
   * read parameters from file
   */
   void HoomdLJShiftedForcePair::readParam(std::istream &in)
   {
      // Read parameters
      readCArray2D<double> (
                  in, "epsilon", epsilon_[0], nAtomType_, nAtomType_);
      readCArray2D<double> (
                  in, "sigma",   sigma_[0], nAtomType_, nAtomType_);
      readCArray2D<double>(
                  in, "cutoff",  cutoff_[0], nAtomType_, nAtomType_);

      // calculate maxPairCutoff and assign parameters
      maxPairCutoff_ = 0.0;
      for (int i = 0; i < nAtomType_; ++i)
         for (int j = 0; j < nAtomType_; ++j) {
            params_[i][j].x = 4.0*epsilon_[i][j]*pow(sigma_[i][j],12.0);
            params_[i][j].y = 4.0*epsilon_[i][j]*pow(sigma_[i][j],6.0);
            if (cutoff_[i][j] > maxPairCutoff_)
               maxPairCutoff_ = cutoff_[i][j];
            cutoffSq_[i][j] = cutoff_[i][j] * cutoff_[i][j];
         }
   }

   void HoomdLJShiftedForcePair::setEpsilon(int i, int j, double epsilon)
   {

      // Preconditions
//      if (!isInitialized_) {
//         UTIL_THROW("Cannot modify epsilon before LJPair is initialized");
//      }
      if (i < 0 || i >= nAtomType_) {
         UTIL_THROW("Invalid atom type index i");
      }
      if (j < 0 || j >= nAtomType_) {
         UTIL_THROW("Invalid atom type index j");
      }

      epsilon_[i][j] = epsilon;
      params_[i][j].x = 4.0*epsilon*pow(sigma_[i][j],12.0);
      params_[i][j].y = 4.0*epsilon*pow(sigma_[i][j],6.0);

      // Symmetrize
      if (j != i) {
         epsilon_[j][i] = epsilon_[i][j];
         params_[j][i] = params_[i][j];
      }
   }

   /* 
   * Get pair interaction strength.
   */
   double HoomdLJShiftedForcePair::epsilon(int i, int j) const
   {
      assert(i >= 0 && i < nAtomType_);
      assert(j >= 0 && j < nAtomType_);
      return epsilon_[i][j];
   }

}

#endif
#endif
