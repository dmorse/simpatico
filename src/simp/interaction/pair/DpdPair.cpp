/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "DpdPair.h"
#ifdef UTIL_MPI
#include <util/mpi/MpiLoader.h>
#endif

#include <iostream>
namespace Simp
{

   using namespace Util;

   /* 
   * Constructor.
   */
   DpdPair::DpdPair() 
    : maxPairCutoff_(0.0),
      nAtomType_(0),
      isInitialized_(false)
   {  setClassName("DpdPair"); }
   
   /* 
   * Copy constructor.
   */
   DpdPair::DpdPair(const DpdPair& other)
    : maxPairCutoff_(other.maxPairCutoff_),
      nAtomType_(other.nAtomType_),
      isInitialized_(other.isInitialized_)
   {
      int i,j;
      for (i = 0; i < nAtomType_; ++i) {
         for (j = 0; j < nAtomType_; ++j) {
            epsilon_[i][j] = other.epsilon_[i][j];
            sigma_[i][j]   = other.sigma_[i][j];
            sigmaSq_[i][j] = other.sigmaSq_[i][j];
            ce_[i][j]      = other.ce_[i][j];
            cf_[i][j]      = other.cf_[i][j];
         } 
      }
   }
   
   /* 
   * Assignment operator.
   */
   DpdPair& DpdPair::operator = (const DpdPair& other)
   {
      maxPairCutoff_ = other.maxPairCutoff_;
      nAtomType_     = other.nAtomType_;
      isInitialized_ = other.isInitialized_;
      int i,j;
      for (i = 0; i < nAtomType_; ++i) {
         for (j = 0; j < nAtomType_; ++j) {
            epsilon_[i][j] = other.epsilon_[i][j];
            sigma_[i][j]   = other.sigma_[i][j];
            sigmaSq_[i][j] = other.sigmaSq_[i][j];
            ce_[i][j]      = other.ce_[i][j];
            cf_[i][j]      = other.cf_[i][j];
         } 
      }
      return *this;
   }
   
   /* 
   * Read potential parameters from file.
   */
   void DpdPair::readParameters(std::istream &in) 
   {
      // Preconditions
      if (nAtomType_ <= 0) {
         UTIL_THROW( "nAtomType must be set before readParam");
      }
   
      // Read parameters
      readCArray2D<double> (in, "epsilon", epsilon_[0], 
                            nAtomType_, nAtomType_, MaxAtomType);
      readCArray2D<double>(in, "sigma", sigma_[0], 
                           nAtomType_, nAtomType_, MaxAtomType);
   
      // Calculate sigmaSq_, ce_, cf_, and maxPairCutoff_
      int i, j;
      maxPairCutoff_ = 0.0;
      for (i = 0; i < nAtomType_; ++i) {
         for (j = 0; j < nAtomType_; ++j) {
            sigmaSq_[i][j]  = sigma_[i][j]*sigma_[i][j];
            if (sigma_[i][j] > maxPairCutoff_) {
               maxPairCutoff_ = sigma_[i][j];
            }
            cf_[i][j] = epsilon_[i][j]/sigmaSq_[i][j];
            ce_[i][j] = 0.5*cf_[i][j];
         }
      }
  
      isInitialized_ = true; 
   }

   /*
   * Load internal state from an archive.
   */
   void DpdPair::loadParameters(Serializable::IArchive &ar)
   {
      // Precondition
      if (nAtomType_ <= 0) {
         UTIL_THROW( "nAtomType must be set before loadParameters");
      }

      // Read parameters
      loadCArray2D<double> (ar, "epsilon", epsilon_[0], 
                            nAtomType_, nAtomType_, MaxAtomType);
      loadCArray2D<double>(ar, "sigma", sigma_[0], 
                            nAtomType_, nAtomType_, MaxAtomType);
      
      #ifdef UTIL_MPI
      MpiLoader<Serializable::IArchive> loader(*this, ar);
      loader.load(sigmaSq_[0], nAtomType_, nAtomType_, MaxAtomType);
      loader.load(cf_[0], nAtomType_, nAtomType_, MaxAtomType);
      loader.load(ce_[0], nAtomType_, nAtomType_, MaxAtomType);
      loader.load(maxPairCutoff_);
      #else
      ar.unpack(sigmaSq_[0], nAtomType_, nAtomType_, MaxAtomType);
      ar.unpack(cf_[0], nAtomType_, nAtomType_, MaxAtomType);
      ar.unpack(ce_[0], nAtomType_, nAtomType_, MaxAtomType);
      ar >> maxPairCutoff_;
      #endif
      isInitialized_ = true;
   }

   /*
   * Save internal state to an archive.
   */
   void DpdPair::save(Serializable::OArchive &ar)
   {
      ar.pack(epsilon_[0], nAtomType_, nAtomType_, MaxAtomType);
      ar.pack(sigma_[0], nAtomType_, nAtomType_, MaxAtomType);
      ar.pack(sigmaSq_[0], nAtomType_, nAtomType_, MaxAtomType);
      ar.pack(cf_[0], nAtomType_, nAtomType_, MaxAtomType);
      ar.pack(ce_[0], nAtomType_, nAtomType_, MaxAtomType);
      ar << maxPairCutoff_;
   }

   /* 
   * Set nAtomType
   */
   void DpdPair::setNAtomType(int nAtomType) 
   {  
      if (nAtomType <= 0) {
         UTIL_THROW("nAtomType <= 0");
      }
      if (nAtomType > MaxAtomType) {
         UTIL_THROW("nAtomType > DpdPair::MaxAtomType");
      }
      nAtomType_ = nAtomType;
   }
    
   /*
   * Reset epsilon_[i][j] after initialization
   */
   void DpdPair::setEpsilon(int i, int j, double epsilon)
   {

      // Preconditions
      if (!isInitialized_) {
         UTIL_THROW("Cannot modify epsilon before DpdPair is initialized");
      }
      if (i < 0 || i >= nAtomType_) {
         UTIL_THROW("Invalid atom type index i");
      }
      if (j < 0 || j >= nAtomType_) {
         UTIL_THROW("Invalid atom type index j");
      }

      epsilon_[i][j] = epsilon;
      cf_[i][j] = epsilon_[i][j]/sigmaSq_[i][j];
      ce_[i][j] = 0.5*cf_[i][j];

      if (j != i) {
         epsilon_[j][i] = epsilon;
         cf_[j][i] = cf_[i][j];
         ce_[j][i] = ce_[i][j];
      }
   } 

   /*
   * Reset sigma_[i][j] after initialization
   */
   void DpdPair::setSigma(int i, int j, double sigma)
   {

      // Preconditions
      if (!isInitialized_) {
         UTIL_THROW("Cannot modify sigma before DpdPair is initialized");
      }
      if (i < 0 || i >= nAtomType_) {
         UTIL_THROW("Invalid atom type index i");
      }
      if (j < 0 || j >= nAtomType_) {
         UTIL_THROW("Invalid atom type index j");
      }

      sigma_[i][j]   = sigma;
      sigmaSq_[i][j] = sigma*sigma;
      cf_[i][j]      = epsilon_[i][j]/sigmaSq_[i][j];
      ce_[i][j]      = 0.5*cf_[i][j];

      // Symmetrize
      if (j != i) {
         sigma_[j][i] = sigma;
         sigmaSq_[j][i] = sigmaSq_[i][j];
         cf_[j][i] = cf_[i][j];
         ce_[j][i] = ce_[i][j];
      }

   } 

   /* 
   * Get pair interaction strength.
   */
   double DpdPair::epsilon(int i, int j) const
   { 
      assert(i >= 0 && i < nAtomType_); 
      assert(j >= 0 && j < nAtomType_); 
      return epsilon_[i][j]; 
   }

   /* 
   * Get pair interaction strength.
   */
   double DpdPair::sigma(int i, int j) const
   {
      assert(i >= 0 && i < nAtomType_); 
      assert(j >= 0 && j < nAtomType_); 
      return sigma_[i][j]; 
   }

   /* 
   * Get maximum of pair cutoff distance, for all atom type pairs.
   */
   double DpdPair::maxPairCutoff() const
   { return maxPairCutoff_; }

   /*
   * Set a potential energy parameter, identified by a string.
   */
   void DpdPair::set(std::string name, int i, int j, double value)
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
   double DpdPair::get(std::string name, int i, int j) const
   {
      double value = 0.0;
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
