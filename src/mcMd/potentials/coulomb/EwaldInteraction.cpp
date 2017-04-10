#ifndef MD_EWALD_INTERACTION_CPP
#define MD_EWALD_INTERACTION_CPP
/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "EwaldInteraction.h"
#ifdef UTIL_MPI
#include <util/mpi/MpiLoader.h>
#endif

#include <iostream>
#include <cstring>
#include <util/math/Constants.h>
#include <cmath>

namespace McMd
{

   using namespace Util;

   /* 
   * Constructor.
   */
   EwaldInteraction::EwaldInteraction() 
    : epsilon_(0.0),
      alpha_(0.0),
      rSpaceCutoff_(0.0),
      kSpaceCutoff_(0.0),
      isInitialized(false)
   { setClassName("EwaldInteraction");}
   
   /* 
   * Copy constructor.
   */
   EwaldInteraction::EwaldInteraction(const EwaldInteraction& other)
    : epsilon_(other.epsilon_),
      alpha_(other.alpha_),
      rSpaceCutoff_(other.rSpaceCutoff_),
      kSpaceCutoff_(other.kSpaceCutoff_),
      isInitialized(other.isInitialized)
   {
      rSpaceCutoffSq_ = other.rSpaceCutoffSq_;
      kSpaceCutoffSq_ = other.kSpaceCutoffSq_;

      /// prefactors for real space energy.
      fourpiepsi_ = 1.0/(epsilon_*4.0*Constants::Pi); 
      twoalpha_ = 2.0*alpha_/sqrt(Constants::Pi);
   }
   
   /* 
   * Assignment operator.
   */
   EwaldInteraction& EwaldInteraction::operator = (const EwaldInteraction& other)
   {
      epsilon_        = other.epsilon_;
      alpha_          = other.alpha_;
      rSpaceCutoff_   = other.rSpaceCutoff_;
      kSpaceCutoff_   = other.kSpaceCutoff_;
      isInitialized  = other.isInitialized;

      rSpaceCutoffSq_ = other.rSpaceCutoffSq_;
      kSpaceCutoffSq_ = other.kSpaceCutoffSq_;

      /// prefactors for real space energy.
      fourpiepsi_ = 1.0/(epsilon_*4.0*Constants::Pi); 
      twoalpha_ = 2.0*alpha_/sqrt(Constants::Pi);
 
      return *this;
   }

   /* 
   * Read potential parameters from file.
   */
   void EwaldInteraction::readParameters(std::istream &in) 
   {
      read<double>(in, "epsilon",      epsilon_);
      read<double>(in, "alpha",        alpha_);
      read<double>(in, "rSpaceCutoff", rSpaceCutoff_);
      read<double>(in, "kSpaceCutoff", kSpaceCutoff_);

      rSpaceCutoffSq_ = rSpaceCutoff_ * rSpaceCutoff_; 
      kSpaceCutoffSq_ = kSpaceCutoff_ * kSpaceCutoff_; 

      /// prefactors for real space energy.
      fourpiepsi_ = 1.0/(epsilon_*4.0*Constants::Pi); 
      twoalpha_ = 2.0*alpha_/sqrt(Constants::Pi);
 
      isInitialized = true;
   }

   /*
   * Load internal state from an archive.
   */
   void EwaldInteraction::loadParameters(Serializable::IArchive &ar)
   {
      // Load all parameters that appear in parameter file
      loadParameter<double>(ar, "epsilon", epsilon_);
      loadParameter<double>(ar, "alpha", alpha_);
      loadParameter<double>(ar, "rSpaceCutoff", rSpaceCutoff_);
      loadParameter<double>(ar, "kSpaceCutoff", kSpaceCutoff_);

      /// Compute prefactors for real space energy and force
      fourpiepsi_ = 1.0/(epsilon_*4.0*Constants::Pi); 
      twoalpha_ = 2.0*alpha_/sqrt(Constants::Pi);
 
      isInitialized = true;
   }

   /*
   * Save internal state to an archive.
   */
   void EwaldInteraction::save(Serializable::OArchive &ar)
   {
      ar << epsilon_;
      ar << alpha_;
      ar << rSpaceCutoff_;
      ar << kSpaceCutoff_;
   }

   /*
   * Modify a parameter, identified by a string.
   */
   void EwaldInteraction::set(std::string name, double value)
   {
      if (name == "epsilon") {
         epsilon_ = value;
      } else
      if (name == "alpha") {
         alpha_ = value;
      } else 
      if (name == "rSpaceCutoff") {
         rSpaceCutoff_ = value;
      } else 
      if (name == "kSpaceCutoff") {
         kSpaceCutoff_ = value;
      } else { 
         UTIL_THROW("Unrecognized parameter name");
      }

      // Recalculate parameter squared.
      rSpaceCutoffSq_ = rSpaceCutoff_ * rSpaceCutoff_;
      kSpaceCutoffSq_ = kSpaceCutoff_ * kSpaceCutoff_;

      /// Compute prefactors for real space energy and force
      fourpiepsi_ = 1.0/(epsilon_*4.0*Constants::Pi); 
      twoalpha_ = 2.0*alpha_/sqrt(Constants::Pi);
   }

   /*
   * Get a parameter value, identified by a string.
   */
   double EwaldInteraction::get(std::string name) const
   {
      double value = 0.0;
      if (name == "epsilon") {
         value = epsilon_;
      } else
      if (name == "alpha") {
         value = alpha_;
      } else
      if (name == "rSpaceCutoff") {
         value = rSpaceCutoff_;
      } else
      if (name == "kSpaceCutoff") {
         value = kSpaceCutoff_;
      } else {
         UTIL_THROW("Unrecognized parameter name");
      }
      return value;
   }

} 
#endif
