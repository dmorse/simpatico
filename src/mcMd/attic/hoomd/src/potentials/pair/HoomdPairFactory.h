#ifndef SIMP_NOPAIR
#ifndef HOOMD_PAIR_FACTORY_H
#define HOOMD_PAIR_FACTORY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <iostream>

#include <mcMd/potentials/pair/PairFactory.h>

#include <hoomd/PotentialPairGPU.h>
#include <hoomd/CellListGPU.h>

#include <modules/hoomd/potentials/pair/HoomdPair.h>

namespace McMd
{

   class System;
   class MdPairPotential;
   class McPairPotential;

   class HoomdPairPotential;

   /**
   * Factory for subclasses MdPairPotential or McPairPotential.
   * 
   * \ingroup Pair_Module
   */
   class HoomdPairFactory : public PairFactory
   {
   
   public:
   
      /**
      * Default constructor.
      */
      HoomdPairFactory();

      // Return simpatico potentials

      /**
      * Return a pointer to a new McPairInteration, if possible.
      */
      McPairPotential* mcFactory(const std::string& subclass, System& system) const;

      /**
      * Return a pointer to a new McPairInteration, if possible.
      */
      MdPairPotential* mdFactory(const std::string& subclass, System& system) const;

      /**
      * Create an MdPairInteraction from a McPairInteration.
      */
      MdPairPotential* mdFactory(McPairPotential& potential) const;

      // Called by HoomdMove

      /**
      * Convert a potential to a HoomdPairPotential pointer
      *
      * \param potential the pair potential
      */
      static HoomdPairPotential *hoomdPairPotentialPtr(McPairPotential&
         potential);
 
      /**
      * Return a ForceCompute
      *
      * \param potential the pair potential the parameters of which are used
      *                  to construct the Hoomd pair potential (must
      *                  be a HoomdPairPotential)
      * \param systemDefinitinoSPtr the Hoomd SystemDefinition
      * \param nListSPtr this will be assigned the Hoomd neighbor list
      * \param skin skin length
      */
      static boost::shared_ptr<ForceCompute>
         hoomdFactory( McPairPotential& potential, System& system,
         boost::shared_ptr<SystemDefinition> systemDefinitionSPtr,
         boost::shared_ptr<NeighborList> &nListSPtr, double skin);

   private:

      /**
      * Internal implementation of hoomdFactory
      */
      template < class hoomd_evaluator,
         cudaError_t gpu_cgpf(const pair_args_t& pair_args,
         const typename hoomd_evaluator::param_type *d_params),
         const char *name>
      static boost::shared_ptr<ForceCompute> hoomdFactoryImpl(
         HoomdPairPotential *hoomdPairPotentialPtr,
         System &sytem,
         boost::shared_ptr<SystemDefinition> systemDefinitionSPtr,
         boost::shared_ptr<NeighborList> &nListSPtr,
         double skin);

   };


   // implementation of the template

   /**
   * return a ForceCompute (implementation)
   */
  template < class hoomd_evaluator,
     cudaError_t gpu_cgpf(const pair_args_t& pair_args,
     const typename hoomd_evaluator::param_type *d_params),
     const char *name>
  boost::shared_ptr<ForceCompute> HoomdPairFactory::hoomdFactoryImpl(
     HoomdPairPotential *hoomdPairPotentialPtr,
     System &system,
     boost::shared_ptr<SystemDefinition> systemDefinitionSPtr,
     boost::shared_ptr<NeighborList> &nListSPtr,
     double skin)
   {
      // NOTE: we are assuming default GPU block sizes (not optimized ones)
      // Create and register cell list
      boost::shared_ptr< CellListGPU > cellListSPtr(
         new CellListGPU(systemDefinitionSPtr));

      HoomdPair < hoomd_evaluator, gpu_cgpf, name >
         *hoomdPairPtr = dynamic_cast< HoomdPair< hoomd_evaluator, gpu_cgpf,
            name> * > (hoomdPairPotentialPtr);

      // Create and register neighbor list
      Scalar defaultRBuff = (Scalar) skin;

      boost::shared_ptr< NeighborListGPUBinned > gpuNList(
         new NeighborListGPUBinned(systemDefinitionSPtr,
         hoomdPairPtr->maxPairCutoff(), defaultRBuff, cellListSPtr));

      // FIXME: need to set optimized values based on GPU capability
      gpuNList->setBlockSize(512);
      gpuNList->setBlockSizeFilter(256);

      nListSPtr = gpuNList;

      // mask bonded pairs ?
      // bonds need to be registered at this stage
      if (system.simulation().maskedPairPolicy() == MaskBonded)
         nListSPtr->addExclusionsFromBonds();

      //nListSPtr->countExclusions();

      // Create and register pair potential
      boost::shared_ptr< PotentialPairGPU< hoomd_evaluator, gpu_cgpf > >
         pairSPtr(new PotentialPairGPU< hoomd_evaluator, gpu_cgpf >(
         systemDefinitionSPtr, nListSPtr,""));

      // Set up energy/force-shift mode
      pairSPtr->setShiftMode( hoomdPairPtr->hoomdShiftMode() );

      // FIXME: need to set optimized values based on GPU capability
      pairSPtr->setBlockSize(192);

      // Set up parameters
      for (int i=0; i < system.simulation().nAtomType(); i++)
         for (int j=0; j< system.simulation().nAtomType(); j++)
         {
            pairSPtr->setParams(i,j,hoomdPairPtr->params(i,j));
            pairSPtr->setRcut(i,j,hoomdPairPtr->cutoff(i,j));
         }

      return pairSPtr;
   }
  
}
#endif
#endif
