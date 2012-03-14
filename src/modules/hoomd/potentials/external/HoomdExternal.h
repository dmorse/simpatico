#ifdef MCMD_EXTERNAL
#ifndef HOOMD_EXTERNAL_H
#define HOOMD_EXTERNAL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/simulation/System.h>
#include <mcMd/simulation/Simulation.h>

#include <util/param/ParamComposite.h>
#include <util/boundary/Boundary.h>
#include <util/global.h>

#include <math.h>

#include <boost/shared_ptr.hpp>

#include <hoomd/ForceCompute.h>
#include <hoomd/System.h>
#include <hoomd/SystemDefinition.h>
#include <hoomd/PotentialExternalGPU.cuh>
#include <hoomd/PotentialExternalGPU.h>
#include <hoomd/HOOMDMath.h>

namespace McMd
{

   using namespace Util;

   /**
   * A potential encapsulating a HOOMD evaluator.
   * Actual implementations have to provide a readParam() method.
   * See HoomdLamellarExternal for an example.
   *
   * \ingroup External_Module
   */
   template < class hoomd_evaluator,
      cudaError_t gpu_cpef(const external_potential_args_t& external_potential_args,
      const typename hoomd_evaluator::param_type *d_params)>
   class HoomdExternal : public ParamComposite
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
      HoomdExternal(const HoomdExternal<hoomd_evaluator, gpu_cpef>& other);

      /**
      * read parameters from file
      */
      virtual void readParam(std::istream &in) = 0;

      /**
      * Set pointer to Boundary.
      *
      * \param boundary Boundary object (used to calculate length along perpendicular direction).
      */
      void setBoundary(Boundary &boundary);
 
      /**
      * Return class name
      */
      virtual std::string className() const = 0;
 
      /* Calculate force divided by distance 
      *
      * \param rsq square of distance between particles
      * \param i   type of particle 1
      * \param j   type of particle 2
      * \return    force divided by distance 
      */
      void getForce(const Vector& position, int type,
                                     Vector& force) const;

      double energy(const Vector& position, int type) const;

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

      /// Are all parameters and pointers initialized?
      bool  isInitialized_;

   private:

      Boundary* boundaryPtr_;
   };
  
   // Implementation

   /*
   * Constructor.
   */
   template < class hoomd_evaluator,
      cudaError_t gpu_cpef(const external_potential_args_t& external_potential_args,
      const typename hoomd_evaluator::param_type *d_params) >
   HoomdExternal< hoomd_evaluator, gpu_cpef >::HoomdExternal()
    : nAtomType_(0),
      isInitialized_(false),
      boundaryPtr_(NULL)
   {
   }    

   /*
   * Copy constructor.
   */
   template < class hoomd_evaluator,
      cudaError_t gpu_cpef(const external_potential_args_t& external_potential_args,
      const typename hoomd_evaluator::param_type *d_params)>
   HoomdExternal<hoomd_evaluator, gpu_cpef>::HoomdExternal(
      const HoomdExternal<hoomd_evaluator, gpu_cpef>& other)
   : nAtomType_(other.nAtomType_),
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
      const typename hoomd_evaluator::param_type *d_params)>
   inline void HoomdExternal< hoomd_evaluator, gpu_cpef >::getForce(
      const Vector& position, int i, Vector& force) const
   {
      Vector lengths = boundaryPtr_->lengths();
      Scalar3 X = make_scalar3(position[0], position[1], position[2]);
      
      hoomd_evaluator eval(X, lengths[0], lengths[1], lengths[2], params_[i]);
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
      const typename hoomd_evaluator::param_type *d_params)>
   inline double HoomdExternal< hoomd_evaluator, gpu_cpef>::energy(
      const Vector& position, int i) const
   {
      Vector lengths = boundaryPtr_->lengths();
      Scalar3 X = make_scalar3(position[0], position[1], position[2]);
      hoomd_evaluator eval(X, lengths[0], lengths[1], lengths[2], params_[i]);
      Scalar3 Force;
      Scalar energy = 0;
      Scalar virial[6];
      eval.evalForceEnergyAndVirial(Force, energy, virial); 
      return energy; 
   }

   /*
   * set the number of atom types
   */
   template < class hoomd_evaluator,
      cudaError_t gpu_cpef(const external_potential_args_t& external_potential_args,
      const typename hoomd_evaluator::param_type *d_params)>
   void HoomdExternal< hoomd_evaluator, gpu_cpef >::setNAtomType(
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
      const typename hoomd_evaluator::param_type *d_params)>
   const typename HoomdExternal<hoomd_evaluator, gpu_cpef>::param_type
    HoomdExternal<hoomd_evaluator, gpu_cpef>::params(int i) const
   {
      return params_[i];
   }

   template < class hoomd_evaluator,
      cudaError_t gpu_cpef(const external_potential_args_t& external_potential_args,
      const typename hoomd_evaluator::param_type *d_params)>
   void HoomdExternal< hoomd_evaluator, gpu_cpef >::setBoundary(Boundary& boundary)
   {
      boundaryPtr_ = &boundary;
   }
}

#endif
#endif
