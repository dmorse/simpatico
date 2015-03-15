#ifdef HOOMD_DEVEL
#ifndef HOOMD_BOND_H
#define HOOMD_BOND_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/mcSimulation/McSystem.h>

#include <util/param/ParamComposite.h>
#include <util/global.h>

#include <math.h>

#include <boost/shared_ptr.hpp>

#include <hoomd/ForceCompute.h>
#include <hoomd/SystemDefinition.h>
#include <hoomd/PotentialBondGPU.h>
#include <hoomd/HOOMDMath.h>

namespace McMd
{

   using namespace Util;

   /*
    * Abstract interface common to all HOOMD bond potentials
    * (used for querying parameters)
    */
   class HoomdBondPotential
   {
   public:
      /**
      * return a ForceCompute for this potential
      *
      * \param systemDefinitionSPtr the HOOMD system 
      * \param system the system
      * \returns a ForceCompute
      */
      virtual boost::shared_ptr<ForceCompute> forceCompute(
         boost::shared_ptr<SystemDefinition> systemDefinitionSPtr,
         McSystem &system) const = 0;

      /**
      * get the internal name of the HOOMD evaluator
      */ 
      virtual std::string evaluatorName() const = 0;
   };

   /**
   * A potential encapsulating a HOOMD evaluator
   *
   * \ingroup Bond_Module
   */
   template < class hoomd_evaluator,
      cudaError_t gpu_cgbf(const bond_args_t& bond_args,
      const typename hoomd_evaluator::param_type *d_params,
      unsigned int *d_flags), const char *name>
   class HoomdBond : public ParamComposite, public HoomdBondPotential
   {
   
   public:
      typedef typename hoomd_evaluator::param_type param_type;
 
      /**
      * Default constructor.
      */
      HoomdBond();

      /**
       * Copy constructor
       */
      HoomdBond(const HoomdBond<hoomd_evaluator, gpu_cgbf, name>& other);

      /**
      * read parameters from file
      */
      void readParam(std::istream &in);
  
      /* Calculate force divided by distance 
      *
      * \param rsq square of distance between particles
      * \param i   type of bond
      * \return    force divided by distance 
      */
      double forceOverR(double rsq, int type) const;

      /* Get interaction energy for a specific bond types
      *
      * \param type type of bond
      * \return bond energy
      */
      double energy(double rsq, int type) const; 

      /**
      * Return class name
      */
      std::string className() const;

      /**  
      * Set nBondType value.
      *
      * \param nAtomType number of atom types.
      */
      void setNBondType(int nBondType);

      /**
       * return a HOOMD ForceCompute
       *
       * \param systemDefinitionSPtr the HOOMD system
       * \param system the system
       */
      virtual boost::shared_ptr<ForceCompute> forceCompute(
         boost::shared_ptr<SystemDefinition> systemDefinitionSPtr,
         McSystem &system) const;

      /**
       * returns internal name of the HOOMD evaluator
       */
      virtual std::string evaluatorName() const;

      /**
      * get the parameters of the HOOMD evaluator
      */
      virtual typename hoomd_evaluator::param_type params(int type) const;

      /**
      * not implemented
      */
      double randomBondLength(Random *random, double beta, int type) const;

   protected:
      /// The encapsulated potential
      boost::shared_ptr<hoomd_evaluator> potential_;  

      /// Number of possible atom types.
      int    nBondType_;

      /// Maximum allowed value for nAtomType (# of particle types)
      static const int MaxNBondType = 2;
 
      /// potential parameters
      param_type params_[MaxNBondType]; 

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
      cudaError_t gpu_cgbf(const bond_args_t& bond_args,
      const typename hoomd_evaluator::param_type *d_params,
      unsigned int *d_flags), const char *name>
   HoomdBond< hoomd_evaluator, gpu_cgbf, name >::HoomdBond()
    : className_(name),
      nBondType_(0),
      isInitialized_(false)
   {
   }    

   /*
   * Copy constructor.
   */
   template < class hoomd_evaluator,
      cudaError_t gpu_cgbf(const bond_args_t& bond_args,
      const typename hoomd_evaluator::param_type *d_params,
      unsigned int *d_flags), const char *name>
   HoomdBond<hoomd_evaluator, gpu_cgbf, name>::HoomdBond(
      const HoomdBond<hoomd_evaluator, gpu_cgbf, name>& other)
   : className_(name),
     nBondType_(other.nBondType_),
     isInitialized_(other.isInitialized_)
   {
      int i,j;
      for (i = 0; i < nBondType_; ++i) {
         params_[i] = other.params_[i];
      }
   }

   /*
   * Get force over distance
   */
   template < class hoomd_evaluator,
      cudaError_t gpu_cgbf(const bond_args_t& bond_args,
      const typename hoomd_evaluator::param_type *d_params,
      unsigned int *d_flags), const char *name>
   inline double HoomdBond< hoomd_evaluator, gpu_cgbf, name >::forceOverR(
      double rsq, int type) const
   {
      hoomd_evaluator eval(rsq, params_[type]);
      Scalar forceOverR =0;
      Scalar bondEnergy = 0;
      eval.evalForceAndEnergy(forceOverR, bondEnergy); 
      return forceOverR; 
   }

   /*
   * Get energy
   */ 
   template < class hoomd_evaluator,
      cudaError_t gpu_cgbf(const bond_args_t& bond_args,
      const typename hoomd_evaluator::param_type *d_params,
      unsigned int *d_flags), const char *name>
   inline double HoomdBond< hoomd_evaluator, gpu_cgbf, name >::energy(
      double rsq, int type) const
   {
      hoomd_evaluator eval(rsq, params_[type]);
      Scalar forceOverR =0;
      Scalar bondEnergy =0;
      eval.evalForceAndEnergy(forceOverR, bondEnergy); 
      return bondEnergy; 
   }

   /*
   * get the class name
   */
   template < class hoomd_evaluator,
      cudaError_t gpu_cgbf(const bond_args_t& bond_args,
      const typename hoomd_evaluator::param_type *d_params,
      unsigned int *d_flags), const char *name>
   inline std::string HoomdBond< hoomd_evaluator, gpu_cgbf, name >::className() const
   {
      return className_;
   }

   /*
   * get internal name of the HOOMD evaluator
   */
   template < class hoomd_evaluator,
      cudaError_t gpu_cgbf(const bond_args_t& bond_args,
      const typename hoomd_evaluator::param_type *d_params,
      unsigned int *d_flags), const char *name>
   inline std::string HoomdBond< hoomd_evaluator, gpu_cgbf, name >
      ::evaluatorName() const
   {
      return hoomd_evaluator::getName();
   }

   /*
   * get the bond parameters
   */
   template < class hoomd_evaluator,
      cudaError_t gpu_cgbf(const bond_args_t& bond_args,
      const typename hoomd_evaluator::param_type *d_params,
      unsigned int *d_flags), const char *name>
   inline typename hoomd_evaluator::param_type HoomdBond<
      hoomd_evaluator, gpu_cgbf, name >::params(int type) const
   {
      return params_[type];
   }

   /*
   * set the number of bond types
   */
   template < class hoomd_evaluator,
      cudaError_t gpu_cgbf(const bond_args_t& bond_args,
      const typename hoomd_evaluator::param_type *d_params,
      unsigned int *d_flags), const char *name>
   void HoomdBond< hoomd_evaluator, gpu_cgbf, name >::setNBondType(
      int nBondType)
   {
      if (nBondType <= 0) {
         UTIL_THROW("nBondType <= 0");
      }
      if (nBondType > MaxNBondType) {
         UTIL_THROW("nAtomType > HoomdBond::MaxAtomType");
      }
      nBondType_ = nBondType;
   }

   /**
   * return a ForceCompute
   */
   template < class hoomd_evaluator,
      cudaError_t gpu_cgbf(const bond_args_t& bond_args,
      const typename hoomd_evaluator::param_type *d_params,
      unsigned int *d_flags), const char *name>
   inline boost::shared_ptr<ForceCompute>
      HoomdBond< hoomd_evaluator, gpu_cgbf, name >::forceCompute(
      boost::shared_ptr<SystemDefinition> systemDefinitionSPtr,
      McSystem &system) const
   {
      // Create and register bond potential
      boost::shared_ptr< PotentialBondGPU< hoomd_evaluator, gpu_cgbf > >
         bondSPtr(new PotentialBondGPU< hoomd_evaluator, gpu_cgbf >(
         systemDefinitionSPtr, ""));

      // FIXME: need to set optimized values based on GPU capability
      bondSPtr->setBlockSize(192);

      // Set up parameters
      for (int i=0; i < nBondType_; i++) {
         bondSPtr->setParams(i,params_[i]);
      }

      return bondSPtr; 
   }

   /* 
   * Read potential parameters from file.
   */
   template < class hoomd_evaluator,
      cudaError_t gpu_cgbf(const bond_args_t& bond_args,
      const typename hoomd_evaluator::param_type *d_params,
      unsigned int *d_flags), const char *name>
   void HoomdBond< hoomd_evaluator, gpu_cgbf, name >::readParam(
      std::istream &in)
   {
      if (nBondType_ <= 0) {
         UTIL_THROW( "nBondType must be set before readParam");
      }

      // Read parameters
      readCArray<param_type> ( in, "params", params_, nBondType_);

      isInitialized_ = true;
   }

   /*
   * Throws exception if called
   */
   template < class hoomd_evaluator,
      cudaError_t gpu_cgbf(const bond_args_t& bond_args,
      const typename hoomd_evaluator::param_type *d_params,
      unsigned int *d_flags), const char *name>
   double HoomdBond< hoomd_evaluator, gpu_cgbf, name >::randomBondLength(
                         Random *random, double beta, int type) const
   {
      UTIL_THROW("Unimplemented function");
   }

}

#endif
#endif
