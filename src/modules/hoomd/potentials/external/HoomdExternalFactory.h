#ifdef INTER_EXTERNAL
#ifndef HOOMD_EXTERNAL_FACTORY_H
#define HOOMD_EXTERNAL_FACTORY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <iostream>

#include <mcMd/potentials/external/ExternalFactory.h>
#include <mcMd/potentials/external/ExternalPotentialImpl.h>

#include <modules/hoomd/potentials/external/HoomdExternal.h>

namespace McMd
{

   class System;

   class HoomdExternalPotential;

   /**
   * Factory for subclasses MdExternalPotential or McExternalPotential.
   * 
   * \ingroup External_Module
   */
   class HoomdExternalFactory : public ExternalFactory
   {
   
   public:
   
      /**
      * Default constructor.
      */
      HoomdExternalFactory(System& system);

      // Return simpatico potentials

      /**
      * Return a pointer to a new ExternalInteration, if possible.
      */
      ExternalPotential* factory(const std::string& subclass) const;

      /**
      * Return a ForceCompute
      *
      * \param potential the external potential the parameters of which are used
      *                  to construct the Hoomd external potential
      * \param system the system
      * \param systemDefinitionSPtr the Hoomd SystemDefinition
      */
      static boost::shared_ptr<ForceCompute>
         hoomdFactory( ExternalPotential& potential, System &system, boost::shared_ptr<SystemDefinition> systemDefinitionSPtr);

   private:

      /**
      * Internal implementation of hoomdFactory
      */
      template < class hoomd_evaluator,
                 cudaError_t gpu_cpef(const external_potential_args_t& external_potential_args,
                                               const typename hoomd_evaluator::param_type *d_params),
                 class T >

      static boost::shared_ptr<ForceCompute> hoomdFactoryImpl(
         ExternalPotential *externalPotentialPtr,
         System &system,
         boost::shared_ptr<SystemDefinition> systemDefinitionSPtr);

      //! reference to system  
      System* systemPtr_;
   };

   /**
   * return a ForceCompute (implementation)
   */
  template < class hoomd_evaluator,
             cudaError_t gpu_cpef(const external_potential_args_t& external_potential_args,
                                               const typename hoomd_evaluator::param_type *d_params),
             class T>
  boost::shared_ptr<ForceCompute> HoomdExternalFactory::hoomdFactoryImpl(
     ExternalPotential *externalPotentialPtr,
     System &system,
     boost::shared_ptr<SystemDefinition> systemDefinitionSPtr)
   {
      ExternalPotentialImpl< T > *evalPtr =
         (dynamic_cast< ExternalPotentialImpl< T > *>( externalPotentialPtr));

      if (! evalPtr)
         UTIL_THROW("Not a Hoomd potential.");

      HoomdExternal < hoomd_evaluator, gpu_cpef >
         *hoomdExternalPtr = dynamic_cast< HoomdExternal<hoomd_evaluator, gpu_cpef> *>(&(evalPtr->interaction()));
      assert(hoomdExternalPtr);
        
      // Create and register external potential
      boost::shared_ptr< PotentialExternalGPU< hoomd_evaluator, gpu_cpef > >
         externalSPtr(new PotentialExternalGPU< hoomd_evaluator, gpu_cpef >(
         systemDefinitionSPtr));


      // FIXME: need to set optimized values based on GPU capability
      //externalSPtr->setBlockSize(512);

      // Set up parameters
      for (int i=0; i < system.simulation().nAtomType(); i++)
      {
         externalSPtr->setParams(i,hoomdExternalPtr->params(i));
      }

      return externalSPtr;
   }
  
}
#endif
#endif
