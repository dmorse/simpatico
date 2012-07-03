#ifndef MCMD_SPECIES_MUTATOR_CPP
#define MCMD_SPECIES_MUTATOR_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "SpeciesMutator.h"
#include <iostream>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   SpeciesMutator::SpeciesMutator()
    : stateIdLabel_("stateId"),
      nState_(0)
   {}

   /*
   * Destructor.
   */
   SpeciesMutator::~SpeciesMutator()
   {} 

   /*
   *
   */
   void SpeciesMutator::allocateSpeciesMutator(int nMolecule, int nState)
   {  
      // Allocate arrays
      stateWeights_.allocate(nState);
      stateOccupancies_.allocate(nState);
      moleculeStateIds_.allocate(nMolecule);
      nState_ = nState;

      // Initialize array elements
      for (int i = 0; i < nState; ++i) {
         stateWeights_[i] = 1.0;
      }
      for (int i = 0; i < nState; ++i) {
         stateOccupancies_[i] = 0;
      }
      for (int i = 0; i < nMolecule; ++i) {
         moleculeStateIds_[i] = NullStateId;
      }

   }

   /*
   * Set the state id for a molecule and update occupancy histogram.
   */
   void SpeciesMutator::setMoleculeStateId(const Molecule& molecule, int stateId)
   {
      if (stateId < 0 && stateId != NullStateId) {
         UTIL_THROW("Invalid state index");
      }
      int molId = molecule.id();
      int oldStateId = moleculeStateIds_[molId];
      if (oldStateId != NullStateId) {
         --stateOccupancies_[oldStateId];
      }
      moleculeStateIds_[molId] = stateId; 
      ++stateOccupancies_[stateId];
   }


   /*
   * Write state index for a molecule to output stream.
   */
   void SpeciesMutator::writeMoleculeState(std::ostream& out, 
		                          const Molecule& molecule) const
   {  
      out << "stateId    "<< moleculeStateId(molecule) << std::endl; 
   }

   /*
   * Read state index for a molecule to output stream, and set molecule state.
   */
   void SpeciesMutator::readMoleculeState(std::istream& in, Molecule& molecule)
   {  
      int stateId;
      in >> stateIdLabel_ >> stateId;
      setMoleculeState(molecule, stateId); 
   }

} 
#endif
