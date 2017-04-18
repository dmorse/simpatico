#ifndef MD_PME_INTERACTION_CPP
#define MD_PME_INTERACTION_CPP
/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "PMEInteraction.h"
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
   PMEInteraction::PMEInteraction() 
    : epsilon_(0.0),
      alpha_(0.0),
      rSpaceCutoff_(0.0),
      xgridSize_(0),
      ygridSize_(0),
      zgridSize_(0),
      isInitialized(false)
   { setClassName("PMEInteraction");}
   
   /* 
   * Copy constructor.
   */
   PMEInteraction::PMEInteraction(const PMEInteraction& other)
    : epsilon_(other.epsilon_),
      alpha_(other.alpha_),
      rSpaceCutoff_(other.rSpaceCutoff_),
      xgridSize_(other.xgridSize_),
      ygridSize_(other.ygridSize_),
      zgridSize_(other.zgridSize_),
      isInitialized(other.isInitialized)
   {
      rSpaceCutoffSq_ = other.rSpaceCutoffSq_;

      /// Prefactors for real space energy.
      ce_ = 1.0/(epsilon_*4.0*Constants::Pi); 
      cf_ = 2.0*alpha_/sqrt(Constants::Pi);
   }
   
   /* 
   * Assignment operator.
   */
   PMEInteraction& PMEInteraction::operator = (const PMEInteraction& other)
   {
      epsilon_ = other.epsilon_;
      alpha_ = other.alpha_;
      xgridSize_ = other.xgridSize_;
      ygridSize_ = other.ygridSize_;
      zgridSize_ = other.zgridSize_;
      rSpaceCutoff_ = other.rSpaceCutoff_;
      isInitialized = other.isInitialized;

      rSpaceCutoffSq_ = other.rSpaceCutoffSq_;

      // Derived constants
      ce_ = 1.0/(epsilon_*4.0*Constants::Pi); 
      cf_ = 2.0*alpha_/sqrt(Constants::Pi);
 
      return *this;
   }

   /* 
   * Read potential parameters from file.
   */
   void PMEInteraction::readParameters(std::istream &in) 
   {
      read<double>(in, "epsilon",      epsilon_);
      read<double>(in, "alpha",        alpha_);
      read<double>(in, "rSpaceCutoff", rSpaceCutoff_);
      read<int>(   in, "xgridSize", xgridSize_);
      read<int>(   in, "ygridSize", ygridSize_);
      read<int>(   in, "zgridSize", zgridSize_);

      // Derived constants
      rSpaceCutoffSq_ = rSpaceCutoff_ * rSpaceCutoff_; 
      ce_ = 1.0/(epsilon_*4.0*Constants::Pi); 
      cf_ = 2.0*alpha_/sqrt(Constants::Pi);
 
      isInitialized = true;
   }

   /*
   * Load internal state from an archive.
   */
   void PMEInteraction::loadParameters(Serializable::IArchive &ar)
   {
      // Load all parameters that appear in parameter file
      loadParameter<double>(ar, "epsilon", epsilon_);
      loadParameter<double>(ar, "alpha", alpha_);
      loadParameter<double>(ar, "rSpaceCutoff", rSpaceCutoff_);
      loadParameter<int>(   ar, "xgridSize", xgridSize_);
      loadParameter<int>(   ar, "ygridSize", ygridSize_);
      loadParameter<int>(   ar, "zgridSize", zgridSize_);

      // Derived constants
      rSpaceCutoffSq_ = rSpaceCutoff_ * rSpaceCutoff_; 
      ce_ = 1.0/(epsilon_*4.0*Constants::Pi); 
      cf_ = 2.0*alpha_/sqrt(Constants::Pi);
 
      isInitialized = true;
   }

   /*
   * Save internal state to an archive.
   */
   void PMEInteraction::save(Serializable::OArchive &ar)
   {
      ar << epsilon_;
      ar << alpha_;
      ar << rSpaceCutoff_;
      ar << xgridSize_;
      ar << ygridSize_;
      ar << zgridSize_;
   }

   /*
   * Modify a parameter, identified by a string.
   */
   void PMEInteraction::set(std::string name, double value)
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
      if (name == "xgridSize") {
         xgridSize_ = value;
      } else  
      if (name == "ygridSize") {
         ygridSize_ = value;
      } else  
      if (name == "zgridSize") {
         zgridSize_ = value;
      } else { 
        UTIL_THROW("Unrecognized parameter name");
      }

      // Recalculate parameter squared.
      rSpaceCutoffSq_ = rSpaceCutoff_ * rSpaceCutoff_;

      /// Compute prefactors for real space energy and force
      ce_ = 1.0/(epsilon_*4.0*Constants::Pi); 
      cf_ = 2.0*alpha_/sqrt(Constants::Pi);
   }

   /*
   * Get a parameter value, identified by a string.
   */
   double PMEInteraction::get(std::string name) const
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
      if (name == "xgridSize") {
         value = xgridSize_;
      } else 
      if (name == "ygridSize") {
         value = ygridSize_;
      } else 
      if (name == "zgridSize") {
         value = zgridSize_;
      } else {
         UTIL_THROW("Unrecognized parameter name");
      }
      return value;
   }

} 
#endif
