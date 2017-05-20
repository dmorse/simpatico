#ifdef SIMP_EXTERNAL
#ifndef HOOMD_EXTERNAL_H
#define HOOMD_EXTERNAL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/simulation/System.h>
#include <mcMd/simulation/Simulation.h>

#include <util/space/Dimension.h>
#include <util/param/ParamComposite.h>
#include <util/boundary/Boundary.h>
#include <util/global.h>

#include <math.h>

#include <boost/shared_ptr.hpp>

#include <hoomd/ForceCompute.h>
#include <hoomd/System.h>
#include <hoomd/SystemDefinition.h>
#include <hoomd/PotentialExternalGPU.h>
#include <hoomd/HOOMDMath.h>
#include <hoomd/BoxDim.h>

#include "HoomdExternalPotential.h"

namespace McMd
{

   using namespace Util;

   /**
   * A potential encapsulating a HOOMD evaluator.
   * Actual implementations have to provide a readParameters() method.
   * See HoomdLamellarOrderingExternal for an example.
   *
   * \ingroup External_Module
   */
   template < class hoomd_evaluator,
      cudaError_t gpu_cpef(const external_potential_args_t& external_potential_args,
      const typename hoomd_evaluator::param_type *d_params),
      const char *name>
   class HoomdExternal : public ParamComposite,
                         public HoomdExternalPotential
   {
   
   public:

      typedef typename hoomd_evaluator::param_type param_type;

      /**
      * Default constructor.
      */
      HoomdExternal();

      /**
      * Copy constructor
      */
      HoomdExternal(const HoomdExternal<hoomd_evaluator, gpu_cpef, name>& other);

      /**
      * read parameters from file
      */
      virtual void readParameters(std::istream &in) = 0;

      /**
      * Set pointer to Boundary.
      *
      * \param boundary Boundary object (used to calculate length along perpendicular direction).
      */
      void setBoundary(Boundary &boundary);
 
      /**
      * Returns internal name of HOOMD evaluator
      */
      virtual std::string hoomdName() const;
 
      /* Calculate force divided by distance 
      *
      * \param rsq square of distance between particles
      * \param i   type of particle 1
      * \param j   type of particle 2
      * \return    force divided by distance 
      */
      void getForce(const Vector& position, int type,
                                     Vector& force) const;

      /* Get interaction energy for a specific atom
      *
      * \param position vector position of atom
      * \param type     type of atom 
      */
      double energy(const Vector& position, int type) const;

      /**
      * Return class name
      */      
      std::string className() const;

      /**  
      * Set nAtomType value.
      *
      * \param nAtomType number of atom types.
      */
      void setNAtomType(int nAtomType);

      /**
      * returns the potential parameters
      * 
      * \param i type of particle
      */
      const param_type params(int i) const;

   protected:

      /// Number of possible atom types.
      int    nAtomType_;

      /// Maximum allowed value for nAtomType (# of particle types)
      static const int MaxAtomType = 2;
 
      /// Potential parameters
      param_type params_[MaxAtomType]; 

   private:

      /// The class name
      std::string className_;

      /// Are all parameters and pointers initialized?
      bool  isInitialized_;

      /// Pointer to system boundary
      Boundary* boundaryPtr_;
   };
  
   // Implementation

   /*
   * Constructor.
   */
   template < class hoomd_evaluator,
      cudaError_t gpu_cpef(const external_potential_args_t& external_potential_args,
      const typename hoomd_evaluator::param_type *d_params),
     const char *name>
   HoomdExternal< hoomd_evaluator, gpu_cpef, name >::HoomdExternal()
    : nAtomType_(0),
      className_(name),
      isInitialized_(false),
      boundaryPtr_(NULL)
   {
   }    

   /*
   * Copy constructor.
   */
   template < class hoomd_evaluator,
      cudaError_t gpu_cpef(const external_potential_args_t& external_potential_args,
      const typename hoomd_evaluator::param_type *d_params),
      const char *name>
   HoomdExternal<hoomd_evaluator, gpu_cpef, name>::HoomdExternal(
      const HoomdExternal<hoomd_evaluator, gpu_cpef, name>& other)
   : nAtomType_(other.nAtomType_),
     className_(name),
     isInitialized_(other.isInitialized_)
   {
      int i;
      for (i = 0; i < nAtomType_; ++i) {
         params_[i] = other.params_[i];
      }
   }


   /*
   * Get force over distance
   */
   template < class hoomd_evaluator,
      cudaError_t gpu_cpef(const external_potential_args_t& external_potential_args,
      const typename hoomd_evaluator::param_type *d_params),
      const char *name>
   inline void HoomdExternal< hoomd_evaluator, gpu_cpef, name >::getForce(
      const Vector& position, int i, Vector& force) const
   {
      Vector lengths = boundaryPtr_->lengths();
      Vector hoomdPosition;
      for (int j=0; j < Dimension; ++j) {
         hoomdPosition[j] = position[j]/lengths[j];
      }
      Scalar3 X = make_scalar3(hoomdPosition[0], hoomdPosition[1], hoomdPosition[2]);
      BoxDim box(lengths[0], lengths[1], lengths[2]);
      hoomd_evaluator eval(X, box, params_[i]);
      Scalar3 Force;
      Scalar energy = 0;
      Scalar virial[6];
      eval.evalForceEnergyAndVirial(Force, energy, virial); 
      force[0] = Force.x; 
      force[1] = Force.y; 
      force[2] = Force.z; 
   }

   /*
   * Get energy
   */ 
   template < class hoomd_evaluator,
      cudaError_t gpu_cpef(const external_potential_args_t& external_potential_args,
      const typename hoomd_evaluator::param_type *d_params),
      const char *name>
   inline double HoomdExternal< hoomd_evaluator, gpu_cpef, name >::energy(
      const Vector& position, int i) const
   {
      Vector lengths = boundaryPtr_->lengths();
      Vector hoomdPosition;
      for (int j=0; j < Dimension; ++j) {
         hoomdPosition[j] = position[j]/lengths[j];
      }
      Scalar3 X = make_scalar3(hoomdPosition[0], hoomdPosition[1], hoomdPosition[2]);
      BoxDim box(lengths[0], lengths[1], lengths[2]);
      hoomd_evaluator eval(X, box, params_[i]);
      Scalar3 Force;
      Scalar energy = 0;
      Scalar virial[6];
      eval.evalForceEnergyAndVirial(Force, energy, virial); 
      return energy; 
   }

   /*
   * get the class name
   */
   template < class hoomd_evaluator,
      cudaError_t gpu_cpef(const external_potential_args_t& external_potential_args,
      const typename hoomd_evaluator::param_type *d_params),
      const char *name>
   std::string HoomdExternal< hoomd_evaluator, gpu_cpef, name >::className() const
   {
      return className_;
   }

   /*
   * get internal name of the HOOMD evaluator
   */
   template < class hoomd_evaluator,
      cudaError_t gpu_cpef(const external_potential_args_t& external_potential_args,
      const typename hoomd_evaluator::param_type *d_params),
      const char *name>
   std::string HoomdExternal< hoomd_evaluator, gpu_cpef, name >
      ::hoomdName() const
   {
      return hoomd_evaluator::getName();
   }

   /*
   * set the number of atom types
   */
   template < class hoomd_evaluator,
      cudaError_t gpu_cpef(const external_potential_args_t& external_potential_args,
      const typename hoomd_evaluator::param_type *d_params),
      const char *name>
   void HoomdExternal< hoomd_evaluator, gpu_cpef, name >::setNAtomType(
      int nAtomType)
   {
      if (nAtomType <= 0) {
         UTIL_THROW("nAtomType <= 0");
      }
      if (nAtomType > MaxAtomType) {
         UTIL_THROW("nAtomType > HoomdPair::MaxAtomType");
      }
      nAtomType_ = nAtomType;
   }

   /**
   * returns the potential parameters
   */
   template < class hoomd_evaluator,
      cudaError_t gpu_cpef(const external_potential_args_t& external_potential_args,
      const typename hoomd_evaluator::param_type *d_params),
      const char *name>
   const typename HoomdExternal<hoomd_evaluator, gpu_cpef, name>::param_type
    HoomdExternal<hoomd_evaluator, gpu_cpef, name>::params(int i) const
   {
      return params_[i];
   }

   template < class hoomd_evaluator,
      cudaError_t gpu_cpef(const external_potential_args_t& external_potential_args,
      const typename hoomd_evaluator::param_type *d_params),
      const char *name>
   void HoomdExternal< hoomd_evaluator, gpu_cpef, name >::setBoundary(Boundary& boundary)
   {
      boundaryPtr_ = &boundary;
   }
}

#endif
#endif
