#ifndef DDMD_SIMULATION_ACCESS_CPP
#define DDMD_SIMULATION_ACCESS_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SimulationAccess.h"
#include "Simulation.h"

namespace DdMd
{

   /*
   * Constructor.
   */
   SimulationAccess::SimulationAccess(Simulation& simulation)
    : 
      atomTypes_(),
      simulationPtr_(&simulation),
      boundaryPtr_(&simulation.boundary_),
      atomStoragePtr_(&simulation.atomStorage_),
      #ifdef INTER_BOND
      bondStoragePtr_(&simulation.bondStorage_),
      #endif
      #ifdef INTER_ANGLE
      angleStoragePtr_(&simulation.angleStorage_),
      #endif
      #ifdef INTER_DIHEDRAL
      dihedralStoragePtr_(&simulation.dihedralStorage_),
      #endif
      pairPotentialPtr_(simulation.pairPotentialPtr_),
      #ifdef INTER_BOND
      bondPotentialPtr_(simulation.bondPotentialPtr_),
      #endif
      #ifdef INTER_ANGLE
      anglePotentialPtr_(simulation.anglePotentialPtr_),
      #endif
      #ifdef INTER_DIHEDRAL
      dihedralPotentialPtr_(simulation.dihedralPotentialPtr_),
      #endif
      #ifdef INTER_EXTERNAL
      externalPotentialPtr_(simulation.externalPotentialPtr_),
      #endif
      energyEnsemblePtr_(simulation.energyEnsemblePtr_),
      boundaryEnsemblePtr_(simulation.boundaryEnsemblePtr_),
      randomPtr_(&simulation.random_),
      domainPtr_(&simulation.domain_),
      exchangerPtr_(&simulation.exchanger_),
      fileMasterPtr_(simulation.fileMasterPtr_),
      nAtomType_(simulation.nAtomType_),
      #ifdef INTER_BOND
      nBondType_(simulation.nBondType_),
      #endif
      #ifdef INTER_ANGLE
      nAngleType_(simulation.nAngleType_),
      #endif
      #ifdef INTER_DIHEDRAL
      nDihedralType_(simulation.nDihedralType_),
      #endif
      #ifdef INTER_EXTERNAL
      hasExternal_(simulation.hasExternal_),
      #endif
      maskedPairPolicy_(simulation.maskedPairPolicy_),
      reverseUpdateFlag_(simulation.reverseUpdateFlag_)
   {
      atomTypes_.associate(simulation.atomTypes_); 
      // RArray atomsTypes_ wraps the same C array as that
      // owned by atomTypes_ DArray in the parent Simulation.
   }

   /*
   * Destructor.
   */
   SimulationAccess::~SimulationAccess()
   {}

}
#endif
