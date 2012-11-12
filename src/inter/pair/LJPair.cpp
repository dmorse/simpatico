#ifndef INTER_LJ_PAIR_CPP
#define INTER_LJ_PAIR_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "LJPair.h"

#include <iostream>
#include <cstring>

namespace Inter
{

   using namespace Util;

   /* 
   * Constructor.
   */
   LJPair::LJPair() 
    : maxPairCutoff_(0.0),
      nAtomType_(0),
      isInitialized_(false)
   {  setClassName("LJPair"); }
   
   /* 
   * Copy constructor.
   */
   LJPair::LJPair(const LJPair& other)
    : maxPairCutoff_(other.maxPairCutoff_),
      nAtomType_(other.nAtomType_),
      isInitialized_(other.isInitialized_)
   {
      int i,j;
      for (i = 0; i < nAtomType_; ++i) {
         for (j = 0; j < nAtomType_; ++j) {
            epsilon_[i][j] = other.epsilon_[i][j];
            sigma_[i][j]   = other.sigma_[i][j];
            cutoff_[i][j]  = other.cutoff_[i][j];
            sigmaSq_[i][j] = other.sigmaSq_[i][j];
            cutoffSq_[i][j]= other.cutoffSq_[i][j];
            ljShift_[i][j] = other.ljShift_[i][j];
         } 
      }
   }
   
   /* 
   * Assignment operator.
   */
   LJPair& LJPair::operator = (const LJPair& other)
   {
      maxPairCutoff_ = other.maxPairCutoff_;
      nAtomType_     = other.nAtomType_;
      isInitialized_ = other.isInitialized_;
      int i,j;
      for (i = 0; i < nAtomType_; ++i) {
         for (j = 0; j < nAtomType_; ++j) {
            epsilon_[i][j] = other.epsilon_[i][j];
            sigma_[i][j]   = other.sigma_[i][j];
            cutoff_[i][j]  = other.cutoff_[i][j];
            sigmaSq_[i][j] = other.sigmaSq_[i][j];
            cutoffSq_[i][j]= other.cutoffSq_[i][j];
            ljShift_[i][j] = other.ljShift_[i][j];
         } 
      }
      return *this;
   }

   /* 
   * Set nAtomType
   */
   void LJPair::setNAtomType(int nAtomType) 
   {  
      if (nAtomType <= 0) {
         UTIL_THROW("nAtomType <= 0");
      }
      if (nAtomType > MaxAtomType) {
         UTIL_THROW("nAtomType > LJPair::MaxAtomType");
      }
      nAtomType_ = nAtomType;
   }

   /*
   * Reset epsilon_[i][j] after initialization
   */
   void LJPair::setEpsilon(int i, int j, double epsilon)
   {

      // Preconditions
      if (!isInitialized_) {
         UTIL_THROW("Cannot modify epsilon before LJPair is initialized");
      }
      if (i < 0 || i >= nAtomType_) {
         UTIL_THROW("Invalid atom type index i");
      }
      if (j < 0 || j >= nAtomType_) {
         UTIL_THROW("Invalid atom type index j");
      }

      epsilon_[i][j] = epsilon;

      // Calculate shift
      double r6i = sigmaSq_[i][j]/cutoffSq_[i][j];
      r6i = r6i*r6i*r6i;
      ljShift_[i][j] = -4.0*epsilon_[i][j]*(r6i*r6i - r6i);

      // Symmetrize
      if (j != i) {
         epsilon_[j][i] = epsilon_[i][j];
         ljShift_[j][i] = ljShift_[i][j];
      }
   } 

   /*
   * Reset sigma_[i][j] after initialization
   */
   void LJPair::setSigma(int i, int j, double sigma)
   {

      // Preconditions
      if (!isInitialized_) {
         UTIL_THROW("Cannot modify sigma before LJPair is initialized");
      }
      if (i < 0 || i >= nAtomType_) {
         UTIL_THROW("Invalid atom type index i");
      }
      if (j < 0 || j >= nAtomType_) {
         UTIL_THROW("Invalid atom type index j");
      }

      sigma_[i][j]   = sigma;
      sigmaSq_[i][j] = sigma*sigma;

      // Recalculate shift
      double r6i = sigmaSq_[i][j]/cutoffSq_[i][j];
      r6i = r6i*r6i*r6i;
      ljShift_[i][j] = -4.0*epsilon_[i][j]*(r6i*r6i - r6i);

      // Symmetrize
      if (j != i) {
         sigma_[j][i]   = sigma_[i][j];
         sigmaSq_[j][i] = sigmaSq_[i][j];
         ljShift_[j][i] = ljShift_[j][i];
      }

   } 

   /* 
   * Read potential parameters from file.
   */
   void LJPair::readParameters(std::istream &in) 
   {
      if (nAtomType_ == 0) {
         UTIL_THROW( "nAtomType must be set before readParam");
      }
   
      // Read parameters
      // //readBegin(in,  "LJPair");
      readCArray2D<double> (
                  in, "epsilon", epsilon_[0], nAtomType_, nAtomType_);
      readCArray2D<double>(
                  in, "sigma", sigma_[0], nAtomType_, nAtomType_);
      readCArray2D<double>(
                  in, "cutoff", cutoff_[0], nAtomType_, nAtomType_);
      //readEnd(in);
   
      // Initialize dependent variables sigmaSq, cutoffSq, and ljShift,
      // and calculate maxPairCutoff_
      double r6i;
      int i, j;
      maxPairCutoff_ = 0.0;
      for (i = 0; i < nAtomType_; ++i) {
         for (j = 0; j < nAtomType_; ++j) {
            sigmaSq_[i][j]  = sigma_[i][j]*sigma_[i][j];
            cutoffSq_[i][j] = cutoff_[i][j]*cutoff_[i][j];
            r6i = sigmaSq_[i][j]/cutoffSq_[i][j];
            r6i = r6i*r6i*r6i;
            ljShift_[i][j] = -4.0*epsilon_[i][j]*(r6i*r6i - r6i);
            if (cutoff_[i][j] > maxPairCutoff_ ) {
               maxPairCutoff_ = cutoff_[i][j];
            }
         } // end for j
      } // end for i

      isInitialized_ = true;
   
   }

   /* 
   * Get maximum of pair cutoff distance, for all atom type pairs.
   */
   double LJPair::maxPairCutoff() const
   { return maxPairCutoff_; }

   /* 
   * Get pair interaction strength.
   */
   double LJPair::epsilon(int i, int j) const
   { 
      assert(i >= 0 && i < nAtomType_); 
      assert(j >= 0 && j < nAtomType_); 
      return epsilon_[i][j]; 
   }

   /* 
   * Get pair interaction strength.
   */
   double LJPair::sigma(int i, int j) const
   {
      assert(i >= 0 && i < nAtomType_); 
      assert(j >= 0 && j < nAtomType_); 
      return sigma_[i][j]; 
   }

   /*
   * Modify a parameter, identified by a string.
   */
   void LJPair::set(std::string name, int i, int j, double value)
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

      sigmaSq_[i][j] = sigma_[i][j]*sigma_[i][j];

      // Recalculate shift
      double r6i = sigmaSq_[i][j]/cutoffSq_[i][j];
      r6i = r6i*r6i*r6i;
      ljShift_[i][j] = -4.0*epsilon_[i][j]*(r6i*r6i - r6i);

      //Symmetrize
      if (j != i) {
         sigmaSq_[j][i] = sigmaSq_[i][j];
         ljShift_[j][i] = ljShift_[j][i];
      }
   }

   /*
   * Get a parameter value, identified by a string.
   */
   double LJPair::get(std::string name, int i, int j) const
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
