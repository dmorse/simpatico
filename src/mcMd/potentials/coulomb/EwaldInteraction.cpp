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
      isInitialized(false)
   { setClassName("EwaldInteraction");}
   
   /* 
   * Copy constructor.
   */
   EwaldInteraction::EwaldInteraction(const EwaldInteraction& other)
    : epsilon_(other.epsilon_),
      alpha_(other.alpha_),
      rSpaceCutoff_(other.rSpaceCutoff_),
      isInitialized(other.isInitialized)
   {
      rSpaceCutoffSq_ = other.rSpaceCutoffSq_;

      /// Derived constants
      ce_ = 1.0/(epsilon_*4.0*Constants::Pi); 
      cf_ = 2.0*alpha_/sqrt(Constants::Pi);
      cg_ = -0.25/(alpha_*alpha_);
   }
   
   /* 
   * Assignment operator.
   */
   EwaldInteraction& EwaldInteraction::operator = (const EwaldInteraction& other)
   {
      epsilon_ = other.epsilon_;
      alpha_ = other.alpha_;
      rSpaceCutoff_ = other.rSpaceCutoff_;
      isInitialized = other.isInitialized;

      rSpaceCutoffSq_ = other.rSpaceCutoffSq_;

      // Derived constants
      ce_ = 1.0/(epsilon_*4.0*Constants::Pi); 
      cf_ = 2.0*alpha_/sqrt(Constants::Pi);
      cg_ = -0.25/(alpha_*alpha_);
 
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

      // Derived constants
      rSpaceCutoffSq_ = rSpaceCutoff_ * rSpaceCutoff_; 
      ce_ = 1.0/(epsilon_*4.0*Constants::Pi); 
      cf_ = 2.0*alpha_/sqrt(Constants::Pi);
      cg_ = -0.25/(alpha_*alpha_);
 
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

      // Derived constants
      rSpaceCutoffSq_ = rSpaceCutoff_ * rSpaceCutoff_; 
      ce_ = 1.0/(epsilon_*4.0*Constants::Pi); 
      cf_ = 2.0*alpha_/sqrt(Constants::Pi);
      cg_ = -0.25/(alpha_*alpha_);
 
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
      } else { 
         UTIL_THROW("Unrecognized parameter name");
      }

      // Recalculate all derived constants
      rSpaceCutoffSq_ = rSpaceCutoff_ * rSpaceCutoff_;
      ce_ = 1.0/(epsilon_*4.0*Constants::Pi); 
      cf_ = 2.0*alpha_/sqrt(Constants::Pi);
      cg_ = -0.25/(alpha_*alpha_);
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
      } else {
         UTIL_THROW("Unrecognized parameter name");
      }
      return value;
   }

} 
#endif
