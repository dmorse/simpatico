#ifndef SIMP_NOPAIR
#ifndef HOOMD_DPD_PAIR_CPP
#define HOOMD_DPD_PAIR_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "HoomdDpdPair.h"


namespace McMd
{
   char classNameHoomdDpd[] = "HoomdDpdPair";

   /**
   * Default constructor.
   */
   HoomdDpdPair::HoomdDpdPair()   
    : HoomdPair< EvaluatorPairDPDThermo, gpu_compute_dpdthermo_forces,
         classNameHoomdDpd >()
   {
   }

   /**
   * Copy constructor
   */
   HoomdDpdPair::HoomdDpdPair(const HoomdDpdPair& other)
    : HoomdPair< EvaluatorPairDPDThermo, gpu_compute_dpdthermo_forces,
         classNameHoomdDpd >(other)
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
   void HoomdDpdPair::readParameters(std::istream &in)
   {

      // Read parameters
      readCArray2D<double> (
                  in, "epsilon", epsilon_[0], nAtomType_, nAtomType_, MaxAtomType);
      readCArray2D<double> (
                  in, "sigma",   sigma_[0], nAtomType_, nAtomType_, MaxAtomType);

      // calculate maxPairCutoff and assign parameters
      maxPairCutoff_ = 0.0;
      for (int i = 0; i < nAtomType_; ++i)
         for (int j = 0; j < nAtomType_; ++j) {
            params_[i][j].x = epsilon_[i][j]; // A parameter
            params_[i][j].y = 0;              // gamma parameter
            cutoff_[i][j] = sigma_[i][j];
            if (cutoff_[i][j] > maxPairCutoff_)
               maxPairCutoff_ = cutoff_[i][j];
            cutoffSq_[i][j] = cutoff_[i][j] * cutoff_[i][j];
         }
   }

   void HoomdDpdPair::setEpsilon(int i, int j, double epsilon)
   {

      // Preconditions
//      if (!isInitialized_) {
//         UTIL_THROW("Cannot modify epsilon before DPDPair is initialized");
//      }
      if (i < 0 || i >= nAtomType_) {
         UTIL_THROW("Invalid atom type index i");
      }
      if (j < 0 || j >= nAtomType_) {
         UTIL_THROW("Invalid atom type index j");
      }

      epsilon_[i][j] = epsilon;
      params_[i][j].x = epsilon;

      // Symmetrize
      if (j != i) {
         epsilon_[j][i] = epsilon_[i][j];
         params_[j][i] = params_[i][j];
      }
   }

   /* 
   * Get pair interaction strength.
   */
   double HoomdDpdPair::epsilon(int i, int j) const
   {
      assert(i >= 0 && i < nAtomType_);
      assert(j >= 0 && j < nAtomType_);
      return epsilon_[i][j];
   }

   /* 
   * Get cutoff distance squared.
   */
   double HoomdDpdPair::cutoffSq(int i, int j) const
   {
      assert(i >= 0 && i < nAtomType_);
      assert(j >= 0 && j < nAtomType_);
      return cutoffSq_[i][j];
   }

   /*
   * Modify a parameter, identified by a string.
   */
   void HoomdDpdPair::set(std::string name, int i, int j, double value)
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
   double HoomdDpdPair::get(std::string name, int i, int j) const
   {
      double value;
      if (name == "epsilon") {
         value = epsilon_[i][j];
      } else
      if (name == "sigma") {
         value = sigma_[i][j];
      } else {
         UTIL_THROW("Unrecognized parameter name");
      }
      return value;
   }

}

#endif
#endif
