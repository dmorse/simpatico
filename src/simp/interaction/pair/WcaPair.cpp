/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "WcaPair.h"

#include <iostream>
#include <cstring>

namespace Simp
{

   using namespace Util;

   /* 
   * Read potential parameters from file.
   */
   void WcaPair::readParameters(std::istream& in) 
   {
      if (nAtomType_ <= 0) {
         UTIL_THROW( "nAtomType must be set before readParam");
      }
   
      // Read parameters epsilon and sigma
      readCArray2D<double>(in, "epsilon", epsilon_[0], 
                           nAtomType_, nAtomType_, MaxAtomType);
      readCArray2D<double>(in, "sigma", sigma_[0], 
                           nAtomType_, nAtomType_, MaxAtomType);
   
      // Initialize all other parameters
      double r6i;
      int i, j;
      maxPairCutoff_ = 0.0;
      for (i = 0; i < nAtomType_; ++i) {
         for (j = 0; j < nAtomType_; ++j) {

            cutoff_[i][j] = sigma_[i][j]*pow(2.0, 1.0/6.0);

            sigmaSq_[i][j] = sigma_[i][j]*sigma_[i][j];
            cutoffSq_[i][j] = cutoff_[i][j]*cutoff_[i][j];
            if (cutoff_[i][j] > maxPairCutoff_ ) {
               maxPairCutoff_ = cutoff_[i][j];
            }

            r6i = sigmaSq_[i][j]/cutoffSq_[i][j];
            r6i = r6i*r6i*r6i;
            ljShift_[i][j] = -4.0*epsilon_[i][j]*(r6i*r6i - r6i);

            eps48_[i][j] = 48.0*epsilon_[i][j];
         } 
      } 
      isInitialized_ = true;
   }

   /*
   * Modify a parameter, identified by a string.
   */
   void WcaPair::set(std::string name, int i, int j, double value)
   {
      if (name == "epsilon") {
         epsilon_[i][j] = value;
         eps48_[i][j] = 48.0*value;
         if (j != i) {
            epsilon_[j][i] = value;
            eps48_[j][i] = eps48_[i][j];
         }
      } else {
         UTIL_THROW("Unrecognized parameter name");
      }

      // Recalculate shift
      double r6i = sigmaSq_[i][j]/cutoffSq_[i][j];
      r6i = r6i*r6i*r6i;
      ljShift_[i][j] = -4.0*epsilon_[i][j]*(r6i*r6i - r6i);
      if (j != i) {
         ljShift_[j][i] = ljShift_[j][i];
      }

   }

} 
