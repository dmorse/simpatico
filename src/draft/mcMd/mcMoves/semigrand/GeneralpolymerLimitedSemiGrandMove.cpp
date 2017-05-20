/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "GeneralpolymerLimitedSemiGrandMove.h"
#include <mcMd/simulation/Simulation.h>
#include <mcMd/mcSimulation/McSystem.h>
#ifndef SIMP_NOPAIR
#include <mcMd/potentials/pair/McPairPotential.h>
#endif
#include <mcMd/species/GeneralpolymerSG.h>
#include <util/global.h>

namespace McMd
{

   using namespace Util;

   /* 
   * Constructor
   */
   GeneralpolymerLimitedSemiGrandMove::GeneralpolymerLimitedSemiGrandMove(McSystem& system) : 
      SystemMove(system),
      speciesId_(-1)
   {  setClassName("GeneralpolymerLimitedSemiGrandMove"); } 
   
   /* 
   * Read parameter speciesId.
   */
   void GeneralpolymerLimitedSemiGrandMove::readParameters(std::istream& in) 
   {
      // Read parameters
      readProbability(in);
      read<int>(in, "speciesId", speciesId_);
      read<int>(in, "UpperLimit", Ulimit_);
      read<int>(in, "LowerLimit", Llimit_);
      // Cast the Species to HomopolymerSG
      speciesPtr_ = dynamic_cast<GeneralpolymerSG*>(&(simulation().species(speciesId_)));
      if (!speciesPtr_) {
         UTIL_THROW("Error: Species must be GeneralpolymerSG");
      }
  
   }
 
   /*
   * Load state from an archive.
   */
   void GeneralpolymerLimitedSemiGrandMove::loadParameters(Serializable::IArchive& ar)
   {  
      McMove::loadParameters(ar);
      loadParameter<int>(ar, "speciesId", speciesId_);

      // Cast the Species to HomopolymerSG
      speciesPtr_ = dynamic_cast<GeneralpolymerSG*>(&(simulation().species(speciesId_)));
      if (!speciesPtr_) {
         UTIL_THROW("Species is not a GeneralpolymerSG");
      }
   }

   /*
   * Save state to an archive.
   */
   void GeneralpolymerLimitedSemiGrandMove::save(Serializable::OArchive& ar)
   {
      McMove::save(ar);
      ar & speciesId_;  
   }

   /* 
   * Generate, attempt and accept or reject a Monte Carlo move.
   */
   bool GeneralpolymerLimitedSemiGrandMove::move() 
   {
      incrementNAttempt();
      Molecule& molecule = system().randomMolecule(speciesId_);
      SpeciesMutator* mutatorPtr = &speciesPtr_->mutator();
      #ifndef SIMP_NOPAIR
      // Calculate pair energy for the chosen molecule
      double oldEnergy = system().pairPotential().moleculeEnergy(molecule);
      #endif

      // Toggle state of the molecule
      int oldStateId = speciesPtr_->mutator().moleculeStateId(molecule);
      int newStateId = (oldStateId == 0) ? 1 : 0;
      speciesPtr_->mutator().setMoleculeState(molecule, newStateId);

      #ifdef SIMP_NOPAIR 

      bool   accept = true;

      #else // ifndef SIMP_NOPAIR
      int newStateTotal = mutatorPtr->stateOccupancy(0);
      // Recalculate pair energy for the molecule
      double newEnergy = system().pairPotential().moleculeEnergy(molecule);

      // Decide whether to accept or reject
      double oldWeight = speciesPtr_->mutator().stateWeight(oldStateId);
      double newWeight = speciesPtr_->mutator().stateWeight(newStateId);
      double ratio  = boltzmann(newEnergy - oldEnergy)*newWeight/oldWeight;
      bool   accept = random().metropolis(ratio);
      #endif

      if (accept && newStateTotal >= Llimit_ && newStateTotal <= Ulimit_) {

         incrementNAccept();

      } else {

         // Revert chosen molecule to original state
         speciesPtr_->mutator().setMoleculeState(molecule, oldStateId);

      }

      return accept;
   }

}
