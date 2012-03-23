#ifndef MCMD_HOMOPOLYMER_SEMI_GRAND_MOVE_CPP
#define MCMD_HOMOPOLYMER_SEMI_GRAND_MOVE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "HomopolymerSemiGrandMove.h"
#include <mcMd/simulation/Simulation.h>
#include <mcMd/mcSimulation/McSystem.h>
#ifndef INTER_NOPAIR
#include <mcMd/potentials/pair/McPairPotential.h>
#endif
#include <mcMd/species/HomopolymerSG.h>
#include <util/global.h>

namespace McMd
{

   using namespace Util;

   /* 
   * Constructor
   */
   HomopolymerSemiGrandMove::HomopolymerSemiGrandMove(McSystem& system) : 
      SystemMove(system),
      speciesId_(-1)
   {} 
   
   /* 
   * Read parameter speciesId.
   */
   void HomopolymerSemiGrandMove::readParam(std::istream& in) 
   {
      // Read parameters
      readProbability(in);
      read<int>(in, "speciesId", speciesId_);

      // Cast the Species to HomopolymerSG
      speciesPtr_ = dynamic_cast<HomopolymerSG*>(&(simulation().species(speciesId_)));
      if (!speciesPtr_) {
         UTIL_THROW("Error: Species must be HomopolymerSG");
      }
  
   }

   /* 
   * Generate, attempt and accept or reject a Monte Carlo move.
   */
   bool HomopolymerSemiGrandMove::move() 
   {
      incrementNAttempt();
      Molecule& molecule = system().randomMolecule(speciesId_);

      #ifndef INTER_NOPAIR
      // Calculate pair energy for the chosen molecule
      double oldEnergy = system().pairPotential().moleculeEnergy(molecule);
      #endif

      // Toggle state of the molecule
      int oldStateId = speciesPtr_->mutator().moleculeStateId(molecule);
      int newStateId = (oldStateId == 0) ? 1 : 0;
      speciesPtr_->mutator().setMoleculeState(molecule, newStateId);

      #ifdef INTER_NOPAIR 

      bool   accept = true;

      #else // ifndef INTER_NOPAIR

      // Recalculate pair energy for the molecule
      double newEnergy = system().pairPotential().moleculeEnergy(molecule);

      // Decide whether to accept or reject
      double oldWeight = speciesPtr_->mutator().stateWeight(oldStateId);
      double newWeight = speciesPtr_->mutator().stateWeight(newStateId);
      double ratio  = boltzmann(newEnergy - oldEnergy)*newWeight/oldWeight;
      bool   accept = random().metropolis(ratio);
   
      #endif

      if (accept) {

         incrementNAccept();

      } else {

         // Revert chosen molecule to original state
         speciesPtr_->mutator().setMoleculeState(molecule, oldStateId);

      }

      return accept;
   }

}
#endif
