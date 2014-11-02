#if defined(MCMD_HOOMD) && defined(HOOMD_DEVEL)
#ifdef MCMD_LINK
#ifndef HOOMD_LINK_BOND_POTENTIAL_H
#define HOOMD_LINK_BOND_POTENTIAL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/


#include <boost/shared_ptr.hpp>
#include <mcMd/potentials/bond/HoomdBond.h>
#include <mcMd/potentials/bond/BondPotentialImpl.h>
#include <mcMd/mcSimulation/McSystem.h>
#include <mcMd/potentials/link/LinkPotentialImpl.h>

namespace McMd
{

   using namespace Util;

   /*
    * Template class for storing a hoomd bond and link potential
    * (used for querying parameters) (their classes are template
    * arguments, and instances are passed to the constructor)
    */
   template < class hoomd_bond_evaluator,
              cudaError_t gpu_bond_cgbf(const bond_args_t& bond_args,
              const typename hoomd_bond_evaluator::param_type *d_params,
              unsigned int *d_flags), const char *bondPotentialName,
              class hoomd_link_evaluator,
              cudaError_t gpu_link_cgbf(const bond_args_t& bond_args,
              const typename hoomd_link_evaluator::param_type *d_params,
              unsigned int *d_flags), const char *linkPotentialName,
              class hoomd_linkbond_evaluator,
              cudaError_t gpu_linkbond_cgbf(const bond_args_t& bond_args,
              const typename hoomd_linkbond_evaluator::param_type *d_params,
              unsigned int *d_flags) >
   class HoomdLinkBondPotentialImpl : public HoomdBondPotential
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
         McSystem &system) const;

      /**
      * get the internal name of the HOOMD evaluator
      */ 
      virtual std::string evaluatorName() const;
   };

   // Implementation

   /*
   * get internal name of the HOOMD evaluator
   */
   template < class hoomd_bond_evaluator,
              cudaError_t gpu_bond_cgbf(const bond_args_t& bond_args,
              const typename hoomd_bond_evaluator::param_type *d_params,
              unsigned int *d_flags), const char *bondPotentialName,
              class hoomd_link_evaluator,
              cudaError_t gpu_link_cgbf(const bond_args_t& bond_args,
              const typename hoomd_link_evaluator::param_type *d_params,
              unsigned int *d_flags), const char *linkPotentialName,
              class hoomd_linkbond_evaluator,
              cudaError_t gpu_linkbond_cgbf(const bond_args_t& bond_args,
              const typename hoomd_linkbond_evaluator::param_type *d_params,
              unsigned int *d_flags) >
   inline std::string HoomdLinkBondPotentialImpl< hoomd_bond_evaluator,
      gpu_bond_cgbf, bondPotentialName, hoomd_link_evaluator, gpu_link_cgbf,
      linkPotentialName, hoomd_linkbond_evaluator, gpu_linkbond_cgbf>
      ::evaluatorName() const
   {
      return hoomd_linkbond_evaluator::getName();
   }

   /**
   * return a ForceCompute
   */
   template < class hoomd_bond_evaluator,
              cudaError_t gpu_bond_cgbf(const bond_args_t& bond_args,
              const typename hoomd_bond_evaluator::param_type *d_params,
              unsigned int *d_flags), const char *bondPotentialName,
              class hoomd_link_evaluator,
              cudaError_t gpu_link_cgbf(const bond_args_t& bond_args,
              const typename hoomd_link_evaluator::param_type *d_params,
              unsigned int *d_flags), const char *linkPotentialName,
              class hoomd_linkbond_evaluator,
              cudaError_t gpu_linkbond_cgbf(const bond_args_t& bond_args,
              const typename hoomd_linkbond_evaluator::param_type *d_params,
              unsigned int *d_flags) >
   inline boost::shared_ptr<ForceCompute>
      HoomdLinkBondPotentialImpl< hoomd_bond_evaluator,
      gpu_bond_cgbf, bondPotentialName, hoomd_link_evaluator, gpu_link_cgbf,
      linkPotentialName, hoomd_linkbond_evaluator, gpu_linkbond_cgbf>
      ::forceCompute( boost::shared_ptr<SystemDefinition> systemDefinitionSPtr,
      McSystem &system) const
   {
      // Create and register bond potential
      boost::shared_ptr< PotentialBondGPU< hoomd_linkbond_evaluator,
         gpu_linkbond_cgbf > >
         linkBondSPtr(new PotentialBondGPU< hoomd_linkbond_evaluator,
         gpu_linkbond_cgbf >(systemDefinitionSPtr, ""));

      // FIXME: need to set optimized values based on GPU capability
      linkBondSPtr->setBlockSize(192);

      // Set up parameters for bonds
      for (int i=0; i < system.simulation().nBondType() ; i++) {
         typename hoomd_linkbond_evaluator::param_type params;
         HoomdBond< hoomd_bond_evaluator, gpu_bond_cgbf,
         bondPotentialName > *bondPtr =
            dynamic_cast< HoomdBond< hoomd_bond_evaluator,
            gpu_bond_cgbf, bondPotentialName > * >(&
            ((dynamic_cast< BondPotentialImpl<HoomdBond< hoomd_bond_evaluator,
            gpu_bond_cgbf, bondPotentialName > > * >(
            &system.bondPotential())->evaluator())));
         assert(bondPtr);
         params.bond_params = bondPtr ->params(i);
         params.type = hoomd_linkbond_evaluator::param_type::is_bond;
         linkBondSPtr->setParams(i,params);
      }

      // Set up parameters for links
      for (int i=0; i < system.simulation().nLinkType(); i++) {
         typename hoomd_linkbond_evaluator::param_type params;
         HoomdBond< hoomd_link_evaluator, gpu_link_cgbf,
         linkPotentialName > *linkPtr =
            dynamic_cast< HoomdBond< hoomd_link_evaluator,
            gpu_link_cgbf, linkPotentialName > * >(&
            ((dynamic_cast< LinkPotentialImpl<HoomdBond< hoomd_link_evaluator,
            gpu_link_cgbf, linkPotentialName > > * >(
            &system.linkPotential()))->evaluator()));
         assert(linkPtr);
         params.link_params = linkPtr->params(i);
         params.type = hoomd_linkbond_evaluator::param_type::is_link;
         // link types are assigned bond types nBondType() <= type <
         // nBondType()+nLinkType()
         linkBondSPtr->setParams(system.simulation().nBondType()+i,params);
      }

      return linkBondSPtr; 
   }
}

#endif
#endif
#endif
