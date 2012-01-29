#ifndef MCMD_NOPAIR
#ifndef PAIR_FACTORY_CPP
#define PAIR_FACTORY_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/potentials/pair/PairFactory.h>

#include <mcMd/simulation/System.h>

// PairPotential interfaces and implementation classes
#include <mcMd/potentials/pair/MdPairPotential.h>
#include <mcMd/potentials/pair/MdPairPotentialImpl.h>
#include <mcMd/potentials/pair/McPairPotential.h>
#include <mcMd/potentials/pair/McPairPotentialImpl.h>

// Pair Potential evaluator classes
#include <mcMd/potentials/pair/LJPair.h>
#include <mcMd/potentials/pair/DpdPair.h>
#include <mcMd/potentials/pair/CompensatedPair.h>

#include <modules/hoomd/potentials/pair/HoomdLJPair.h>
#include <modules/hoomd/potentials/pair/HoomdLJShiftedForcePair.h>
#include <modules/hoomd/potentials/pair/HoomdDpdPair.h>

#include <mcMd/potentials/bond/FeneBond.h>

#include <hoomd/HOOMDMath.h>
#include <hoomd/EvaluatorPairLJ.h>
#include <hoomd/EvaluatorPairDPDThermo.h>
#include <hoomd/AllDriverPotentialPairGPU.cuh>

#include "HoomdPairFactory.h"

namespace McMd
{

   /**
   * Default constructor.
   */
   HoomdPairFactory::HoomdPairFactory()
   {}

   /*
   * Return a pointer to a new McPairInteration, if possible.
   */
   McPairPotential* 
   HoomdPairFactory::mcFactory(const std::string& name, System& system) const
   {
      McPairPotential* ptr = 0;
      if (name == classNameHoomdLJ) {
         ptr = new McPairPotentialImpl< HoomdLJPair> (system);
      } else
      if (name == classNameHoomdLJShiftedForce) {
         ptr = new McPairPotentialImpl< HoomdLJShiftedForcePair > (system);
      } else
      if (name == classNameHoomdDpd) {
         ptr = new McPairPotentialImpl< HoomdDpdPair >(system);
      }
      return ptr;
   }

   /*
   * Return a pointer to a new MdPairPotential, if possible.
   */
   MdPairPotential* 
   HoomdPairFactory::mdFactory(const std::string& name, System& system) const
   {
      MdPairPotential* ptr = 0;
      
      if (name == classNameHoomdLJ) {
         ptr = new MdPairPotentialImpl< HoomdLJPair >(system);
      } else if (name == classNameHoomdLJShiftedForce) {
         ptr = new MdPairPotentialImpl< HoomdLJShiftedForcePair >(system);
      } else if (name == classNameHoomdDpd) {
         ptr = new MdPairPotentialImpl< HoomdDpdPair >(system);
      }

      return ptr;
   }

   /*
   * Convert an McPairPotential to a MdPairPotential, if possible.
   */
   MdPairPotential* 
   HoomdPairFactory::mdFactory(McPairPotential& potential) const
   {
      std::string name = potential.evaluatorClassName();
      MdPairPotential* ptr = 0;
      if (name == classNameHoomdLJ) {
         McPairPotentialImpl< HoomdLJPair > *mcPtr
            = dynamic_cast< McPairPotentialImpl< HoomdLJPair > * >(&potential);
         ptr = new MdPairPotentialImpl< HoomdLJPair > (*mcPtr);
      } else if (name == classNameHoomdLJShiftedForce) {
         McPairPotentialImpl< HoomdLJShiftedForcePair > *mcPtr
            = dynamic_cast< McPairPotentialImpl< HoomdLJShiftedForcePair > * >
            (&potential);
         ptr = new MdPairPotentialImpl< HoomdLJShiftedForcePair > (*mcPtr);
      } else if (name == classNameHoomdDpd) {
         McPairPotentialImpl< HoomdDpdPair > *mcPtr
            = dynamic_cast< McPairPotentialImpl< HoomdDpdPair > * >(&potential);
         ptr = new MdPairPotentialImpl< HoomdDpdPair >(*mcPtr);
      }

      return ptr;
   }

   /** 
   * return a HoomdPairPotential
   */
   HoomdPairPotential *HoomdPairFactory::hoomdPairPotentialPtr(McPairPotential&
      potential)
   {
      std::string name = potential.evaluatorClassName();
      HoomdPairPotential* ptr = 0;

      if (name == classNameHoomdLJ) {
         ptr = dynamic_cast< HoomdPairPotential *>(
            dynamic_cast< HoomdLJPair * >(
            &(dynamic_cast< McPairPotentialImpl< HoomdLJPair > * >
            (&potential))->evaluator()));
      } else if (name == classNameHoomdLJShiftedForce) {
         ptr = dynamic_cast< HoomdPairPotential *>(
            dynamic_cast< HoomdLJShiftedForcePair * >(
            &(dynamic_cast< McPairPotentialImpl< HoomdLJShiftedForcePair > * >
            (&potential))->evaluator()));
      } else if (name == classNameHoomdDpd) {
         ptr = dynamic_cast<HoomdPairPotential *>(
            dynamic_cast< HoomdDpdPair* >(
            &(dynamic_cast< McPairPotentialImpl< HoomdDpdPair > * >
            (&potential))->evaluator()));
      }
      return ptr;
   }


   /**
   * return a ForceCompute (interface)
   */
   boost::shared_ptr<ForceCompute> HoomdPairFactory::hoomdFactory(
      McPairPotential& potential,
      System &system,
      boost::shared_ptr<SystemDefinition> systemDefinitionSPtr,
      boost::shared_ptr<NeighborList> &nListSPtr, double skin)
   {
      // first get a pointer to a HoomdPairPotential
      HoomdPairPotential* ptr = hoomdPairPotentialPtr(potential);

      std::string className = potential.evaluatorClassName();
      boost::shared_ptr<ForceCompute> pairSPtr;

      if (className == classNameHoomdLJ) {
         pairSPtr = hoomdFactoryImpl<EvaluatorPairLJ, gpu_compute_ljtemp_forces,
            classNameHoomdLJ>(ptr,system,systemDefinitionSPtr, nListSPtr, skin);
      } else if (className == classNameHoomdLJShiftedForce) {
         pairSPtr = hoomdFactoryImpl<EvaluatorPairLJ, gpu_compute_ljtemp_forces,
            classNameHoomdLJShiftedForce>(ptr,system,systemDefinitionSPtr,
            nListSPtr, skin);
      } else if (className == classNameHoomdDpd) {
         pairSPtr = hoomdFactoryImpl<EvaluatorPairDPDThermo,
            gpu_compute_dpdthermo_forces, classNameHoomdDpd>(ptr,system,
            systemDefinitionSPtr, nListSPtr, skin);
      } else
         UTIL_THROW("Unsupported Hoomd potential." );

      return pairSPtr; 
   }


}
#endif
#endif
