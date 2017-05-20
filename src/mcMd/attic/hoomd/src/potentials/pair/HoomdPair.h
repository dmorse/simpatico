#ifndef SIMP_NOPAIR
#ifndef HOOMD_PAIR_H
#define HOOMD_PAIR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/simulation/System.h>
#include <mcMd/simulation/Simulation.h>

#include <util/param/ParamComposite.h>
#include <util/global.h>

#include <math.h>

#include <boost/shared_ptr.hpp>

#include <hoomd/ForceCompute.h>
#include <hoomd/CellListGPU.h>
#include <hoomd/NeighborListGPUBinned.h>
#include <hoomd/System.h>
#include <hoomd/SystemDefinition.h>
#include <hoomd/PotentialPairGPU.h>
#include <hoomd/HOOMDMath.h>

#include "HoomdPairPotential.h"

namespace McMd
{

   using namespace Util;

   /**
   * A potential encapsulating a HOOMD evaluator.
   * Actual implementations have to provide a readParameters() method.
   * See HoomdLJPair for an example.
   *
   * \ingroup Pair_Module
   */
   template < class hoomd_evaluator,
      cudaError_t gpu_cgpf(const pair_args_t& pair_args,
      const typename hoomd_evaluator::param_type *d_params),
      const char *name>
   class HoomdPair : public ParamComposite, public HoomdPairPotential
   {
   
   public:

      typedef typename hoomd_evaluator::param_type param_type;
 
      /**
      * Default constructor.
      */
      HoomdPair();

      /**
      * Copy constructor
      */
      HoomdPair(const HoomdPair<hoomd_evaluator, gpu_cgpf, name>& other);

      /**
      * read parameters from file
      */
      void readParameters(std::istream &in);
  
      /* Calculate force divided by distance 
      *
      * \param rsq square of distance between particles
      * \param i   type of particle 1
      * \param j   type of particle 2
      * \return    force divided by distance 
      */
      double forceOverR(double rsq, int i, int j) const;

      /* Get interaction energy for a specific pair of Atom types
      *
      * \param i   type of Atom 1
      * \param j   type of Atom 2
      * \return    epsilon_[i][j]
      */
      double energy(double rsq, int i, int j) const; 

      /**
      * Return class name
      */
      std::string className() const;

      /**
      * Get maximum of pair cutoff distance, for all atom type pairs.
      */
      double maxPairCutoff() const;

      /**  
      * Set nAtomType value.
      *
      * \param nAtomType number of atom types.
      */
      void setNAtomType(int nAtomType);

      /**
      * returns internal name of the HOOMD evaluator
      */
      virtual std::string hoomdName() const;

      /**
      * returns the potential parameters
      * 
      * \param i type of first particle
      * \param j type of second particle
      */
      const param_type params(int i, int j) const;

      /**
      * returns the cutoff
      *
      * \param i type of first particle
      * \param j type of second particle
      */
      const double cutoff(int i, int j) const;

      /**
      * Get the Hoomd shift mode for this potential.
      */
      virtual typename PotentialPairGPU< hoomd_evaluator, gpu_cgpf >::
         energyShiftMode hoomdShiftMode() const = 0;

   protected:
      /**
      * Maximum pair potential cutoff radius, for all monomer type pairs.
      *
      * Used in construction of a cell list or Verlet pair list.
      */
      double maxPairCutoff_;

      /// Number of possible atom types.
      int    nAtomType_;

      /// Maximum allowed value for nAtomType (# of particle types)
      static const int MaxAtomType = 2;
 
      /// Cutoff distance.
      double cutoff_[MaxAtomType][MaxAtomType];

      /// Squared cutoff distance
      double cutoffSq_[MaxAtomType][MaxAtomType];   

      /// Potential parameters
      param_type params_[MaxAtomType][MaxAtomType]; 

   private:

      /// The class name
      std::string className_;

      /// Are all parameters and pointers initialized?
      bool  isInitialized_;

   };
  
   // Implementation

   /*
   * Constructor.
   */
   template < class hoomd_evaluator,
      cudaError_t gpu_cgpf(const pair_args_t& pair_args,
      const typename hoomd_evaluator::param_type *d_params),
      const char *name>
   HoomdPair< hoomd_evaluator, gpu_cgpf, name >::HoomdPair()
    : maxPairCutoff_(0.0),
      nAtomType_(0),
      className_(name),
      isInitialized_(false)
   {
   }    

   /*
   * Copy constructor.
   */
   template < class hoomd_evaluator,
      cudaError_t gpu_cgpf(const pair_args_t& pair_args,
      const typename hoomd_evaluator::param_type *d_params),
      const char *name>
   HoomdPair<hoomd_evaluator, gpu_cgpf, name>::HoomdPair(
      const HoomdPair<hoomd_evaluator, gpu_cgpf, name>& other)
   : maxPairCutoff_(other.maxPairCutoff_),
     nAtomType_(other.nAtomType_),
     className_(name),
     isInitialized_(other.isInitialized_)
   {
      int i,j;
      for (i = 0; i < nAtomType_; ++i) {
         for (j = 0; j < nAtomType_; ++j) {
            cutoff_[i][j] = other.cutoff_[i][j];
            cutoffSq_[i][j] = other.cutoffSq_[i][j];
            params_[i][j] = other.params_[i][j];
         }
      }
   }

   /*
   * Get force over distance
   */
   template < class hoomd_evaluator,
      cudaError_t gpu_cgpf(const pair_args_t& pair_args,
      const typename hoomd_evaluator::param_type *d_params),
      const char *name>
   inline double HoomdPair< hoomd_evaluator, gpu_cgpf, name >::forceOverR(
      double rsq, int i, int j) const
   {
      hoomd_evaluator eval(rsq, cutoffSq_[i][j], params_[i][j]);
      Scalar forceOverR =0;
      Scalar pairEnergy = 0;
      bool energyShift = true;
      eval.evalForceAndEnergy(forceOverR, pairEnergy, energyShift); 
      return forceOverR; 
   }

   /*
   * Get energy
   */ 
   template < class hoomd_evaluator,
      cudaError_t gpu_cgpf(const pair_args_t& pair_args,
      const typename hoomd_evaluator::param_type *d_params),
      const char *name>
   inline double HoomdPair< hoomd_evaluator, gpu_cgpf, name >::energy(
      double rsq, int i, int j) const
   {
      hoomd_evaluator eval(rsq, cutoffSq_[i][j], params_[i][j]);
      Scalar forceOverR =0;
      Scalar pairEnergy =0;
      bool energyShift = true;
      eval.evalForceAndEnergy(forceOverR, pairEnergy, energyShift); 
      return pairEnergy; 
   }

   /*
   * get the class name
   */
   template < class hoomd_evaluator,
      cudaError_t gpu_cgpf(const pair_args_t& pair_args,
      const typename hoomd_evaluator::param_type *d_params),
      const char *name>
   std::string HoomdPair< hoomd_evaluator, gpu_cgpf, name >::className() const
   {
      return className_;
   }

   /*
   * get internal name of the HOOMD evaluator
   */
   template < class hoomd_evaluator,
      cudaError_t gpu_cgpf(const pair_args_t& pair_args,
      const typename hoomd_evaluator::param_type *d_params),
      const char *name>
   std::string HoomdPair< hoomd_evaluator, gpu_cgpf, name >
      ::hoomdName() const
   {
      return hoomd_evaluator::getName();
   }

   template < class hoomd_evaluator,
      cudaError_t gpu_cgpf(const pair_args_t& pair_args,
      const typename hoomd_evaluator::param_type *d_params),
      const char *name>
   double  HoomdPair< hoomd_evaluator, gpu_cgpf, name >::maxPairCutoff() const
   { return maxPairCutoff_; }

   /*
   * set the number of atom types
   */
   template < class hoomd_evaluator,
      cudaError_t gpu_cgpf(const pair_args_t& pair_args,
      const typename hoomd_evaluator::param_type *d_params),
      const char *name>
   void HoomdPair< hoomd_evaluator, gpu_cgpf, name >::setNAtomType(
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

   /* 
   * Read potential parameters from file.
   */
   template < class hoomd_evaluator,
      cudaError_t gpu_cgpf(const pair_args_t& pair_args,
      const typename hoomd_evaluator::param_type *d_params),
      const char *name>
   void HoomdPair< hoomd_evaluator, gpu_cgpf, name >::readParameters(
      std::istream &in)
   {
      UTIL_THROW("readParameters() not implemented.");
   }

   /**
   * returns the potential parameters
   */
   template < class hoomd_evaluator,
      cudaError_t gpu_cgpf(const pair_args_t& pair_args,
      const typename hoomd_evaluator::param_type *d_params),
      const char *name>
   const typename HoomdPair<hoomd_evaluator, gpu_cgpf, name>::param_type
    HoomdPair<hoomd_evaluator, gpu_cgpf, name>::params(int i, int j) const
   {
      return params_[i][j];
   }

   /**
   * returns the cutoff
   */
   template < class hoomd_evaluator,
      cudaError_t gpu_cgpf(const pair_args_t& pair_args,
      const typename hoomd_evaluator::param_type *d_params),
      const char *name>
   const double HoomdPair<hoomd_evaluator, gpu_cgpf, name>::
      cutoff(int i, int j) const
   {
      return cutoff_[i][j];
   }

}

#endif
#endif
