#ifndef END_SWAP_MOVE_CPP
#define MCMD_END_SWAP_MOVE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "EndSwapMove.h"
#include <mcMd/simulation/Simulation.h>
#include <mcMd/mcSimulation/McSystem.h>
#include <mcMd/mcSimulation/mc_potentials.h>
#include <util/boundary/Boundary.h>
#include <mcMd/species/Species.h>
#include <mcMd/species/Linear.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <util/global.h>

namespace McMd
{

   using namespace Util;

   /* 
   * Constructor
   */
   EndSwapMove::EndSwapMove(McSystem& system) : 
      SystemMove(system),
      speciesId_(-1)
   {} 
   
   /* 
   * Read parameter speciesId.
   */
   void EndSwapMove::readParam(std::istream& in) 
   {
      readProbability(in);
      read<int>(in, "speciesId", speciesId_);

      Species* speciesPtr = &(simulation().species(speciesId_));
      int nAtom = speciesPtr->nAtom();

      // Preconditions: species must be immutable and a subclass of Linear 
      if (speciesPtr->isMutable()) {
         UTIL_THROW("EndSwapMove applied to mutable species");
      }
      Linear* linearPtr = dynamic_cast<Linear*>(speciesPtr);
      if (linearPtr == 0) {
         UTIL_THROW("EndSwapMove applied to species that is not Linear");
      }
  
      // Allocate memory 
      atomTypeIds_.allocate(nAtom);
      positions_.allocate(nAtom);

      // Set array of atom type ids
      for (int i = 0; i < nAtom; ++i) {
           atomTypeIds_[i] = speciesPtr->atomTypeId(i); 
      }

   }

   /* 
   * Generate, attempt and accept or reject a Monte Carlo move.
   */
   bool EndSwapMove::move() 
   {
      double newEnergy, oldEnergy;
      Molecule* molPtr;
      Atom* atomPtr;
      int i, nAtom;

      incrementNAttempt();

      molPtr = &(system().randomMolecule(speciesId_));
      nAtom  = molPtr->nAtom();

      // Calculate old molecule energy = pair + external 
      oldEnergy = 0.0;
      #ifndef INTER_NOPAIR 
      oldEnergy += system().pairPotential().moleculeEnergy(*molPtr);
      #endif
      #ifdef INTER_EXTERNAL
      for (i = 0; i < nAtom; ++i) {
         atomPtr = &molPtr->atom(i);
         oldEnergy += system().externalPotential().atomEnergy(*atomPtr);
      }
      #endif

      // Reverse sequence of atom types
      for (i = 0; i < nAtom; ++i) {
         assert(molPtr->atom(i).typeId() == atomTypeIds_[i]);
         molPtr->atom(i).setTypeId(atomTypeIds_[nAtom - 1 - i]);   
      }

      // Calculate new energy (with reversed sequence of atom types).
      newEnergy = 0.0;
      #ifndef INTER_NOPAIR 
      newEnergy += system().pairPotential().moleculeEnergy(*molPtr);
      #endif
      #ifdef INTER_EXTERNAL
      for (i = 0; i < nAtom; ++i) {
         molPtr->atom(i).setTypeId(atomTypeIds_[nAtom - 1 - i]);   
         atomPtr = &molPtr->atom(i);
         newEnergy += system().externalPotential().atomEnergy(*atomPtr);
      }
      #endif

      // Restore original sequence of atom type Ids
      for (i = 0; i < nAtom; ++i) {
         molPtr->atom(i).setTypeId(atomTypeIds_[i]);   
      }
   
      // Decide whether to accept or reject
      bool accept = random().metropolis(boltzmann(newEnergy-oldEnergy));
 
      if (accept) {

         // Store all atomic positions
         for (i = 0; i < nAtom; ++i) {
            positions_[i] = molPtr->atom(i).position();
         }

         // Reverse sequence of atomic positions 
         for (i = 0; i < nAtom; ++i) {
            atomPtr = &molPtr->atom(i);
            atomPtr->position() = positions_[nAtom - 1 - i];
            #ifndef INTER_NOPAIR
            system().pairPotential().updateAtomCell(*atomPtr);
            #endif
         }

         incrementNAccept();

      }
      return accept;
   
   }

}
#endif
