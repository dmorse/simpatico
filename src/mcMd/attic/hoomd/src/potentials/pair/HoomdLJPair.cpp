#ifndef SIMP_NOPAIR
#ifndef HOOMD_LJ_PAIR_CPP
#define HOOMD_LJ_PAIR_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "HoomdLJPair.h"


namespace McMd
{
   char classNameHoomdLJ[] = "HoomdLJPair";

   /**
   * Default constructor.
   */
   HoomdLJPair::HoomdLJPair()   
    : HoomdPair< EvaluatorPairLJ, gpu_compute_ljtemp_forces,
         classNameHoomdLJ >()
   {
   }

   /**
   * Copy constructor
   */
   HoomdLJPair::HoomdLJPair(const HoomdLJPair& other)
    : HoomdPair< EvaluatorPairLJ, gpu_compute_ljtemp_forces,
         classNameHoomdLJ >(other)
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
   void HoomdLJPair::readParameters(std::istream &in)
   {
      // Read parameters
      readCArray2D<double> (
                  in, "epsilon", epsilon_[0], nAtomType_, nAtomType_, MaxAtomType);
      readCArray2D<double> (
                  in, "sigma",   sigma_[0], nAtomType_, nAtomType_, MaxAtomType);
      readCArray2D<double>(
                  in, "cutoff",  cutoff_[0], nAtomType_, nAtomType_, MaxAtomType);

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

   void HoomdLJPair::setEpsilon(int i, int j, double epsilon)
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
   double HoomdLJPair::epsilon(int i, int j) const
   {
      assert(i >= 0 && i < nAtomType_);
      assert(j >= 0 && j < nAtomType_);
      return epsilon_[i][j];
   }

   /* 
   * Get cutoff distance squared.
   */
   double HoomdLJPair::cutoffSq(int i, int j) const
   {
      assert(i >= 0 && i < nAtomType_);
      assert(j >= 0 && j < nAtomType_);
      return cutoffSq_[i][j];
   }

   /*
   * Modify a parameter, identified by a string.
   */
   void HoomdLJPair::set(std::string name, int i, int j, double value)
   {
      if (name == "epsilon") {
         epsilon_[i][j] = value;
         epsilon_[j][i] = value;
      } else
      if (name == "sigma") {
         sigma_[i][j] = value;
         sigma_[j][i] = value;
      } else {
         UTIL_THROW("Unrecognized parameter name");
      }
   }

   /*
   * Get a parameter value, identified by a string.
   */
   double HoomdLJPair::get(std::string name, int i, int j) const
   {
      double value;
      if (name == "epsilon") {
         value = epsilon_[i][j];
      } else
      if (name == "sigma") {
         value = sigma_[i][j];
      } else
      if (name == "cutoff") {
         value = cutoff_[i][j];
      } else {
         UTIL_THROW("Unrecognized parameter name");
      }
      return value;
   }

}

#endif
#endif
