/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "EndSwapMove.h"
#include <mcMd/simulation/Simulation.h>
#include <mcMd/mcSimulation/McSystem.h>
#include <mcMd/mcSimulation/mc_potentials.h>
#include <util/boundary/Boundary.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>

#include <simp/species/Species.h>
#include <simp/species/Linear.h>

#include <util/global.h>

namespace McMd
{

   using namespace Util;
   using namespace Simp;

   /* 
   * Constructor
   */
   EndSwapMove::EndSwapMove(McSystem& system) : 
      SystemMove(system),
      speciesId_(-1)
   {
      /* 
      * Preconditions:  
      * For now, Usage with angles and dihedrals is prohibited. This
      * move would work fine with homogeneous angles and dihedrals,
      * and could be modified to work with heterogeneous potentials,
      * but the required checks or modifications are not implemented.
      */

      #ifdef SIMP_ANGLE
      if (system.hasAnglePotential()) {
         UTIL_THROW("CfbEndBase unusable with heterogeneous dihedrals");
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (system.hasDihedralPotential()) {
         UTIL_THROW("CfbEndBase unusable with heterogeneous dihedrals");
      }
      #endif

      setClassName("EndSwapMove"); 
   }
   
   /* 
   * Read parameter speciesId.
   */
   void EndSwapMove::readParameters(std::istream& in) 
   {
      readProbability(in);
      read<int>(in, "speciesId", speciesId_);

      Species* speciesPtr = &(simulation().species(speciesId_));
      int nAtom = speciesPtr->nAtom();

      // Preconditions
      if (speciesPtr->isMutable()) {
         UTIL_THROW("EndSwapMove on mutable Species");
      }
      Linear* linearPtr = dynamic_cast<Linear*>(speciesPtr);
      if (linearPtr == 0) {
         UTIL_THROW("EndSwapMove on Species that is not a Linear");
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
   * Load from archive.
   */
   void EndSwapMove::loadParameters(Serializable::IArchive& ar) 
   {
      McMove::loadParameters(ar);
      loadParameter<int>(ar, "speciesId", speciesId_);
      ar & atomTypeIds_;

      // Validate
      Species* speciesPtr = &(simulation().species(speciesId_));
      int nAtom = speciesPtr->nAtom();
      if (speciesPtr->isMutable()) {
         UTIL_THROW("EndSwapMove applied to mutable species");
      }
      Linear* linearPtr = dynamic_cast<Linear*>(speciesPtr);
      if (linearPtr == 0) {
         UTIL_THROW("EndSwapMove applied to species that is not Linear");
      }
      if (nAtom != atomTypeIds_.capacity()) {
         UTIL_THROW("Inconsistent capacity for atomTypeIds array");
      }
  
      positions_.allocate(nAtom);
   }

   /* 
   * Save to archive.
   */
   void EndSwapMove::save(Serializable::OArchive& ar) 
   {
      McMove::save(ar);
      ar & speciesId_;
      ar & atomTypeIds_;
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
      #ifndef SIMP_NOPAIR 
      oldEnergy += system().pairPotential().moleculeEnergy(*molPtr);
      #endif
      #ifdef SIMP_EXTERNAL
      if (system().hasExternalPotential()) {
         for (i = 0; i < nAtom; ++i) {
            atomPtr = &molPtr->atom(i);
            oldEnergy += system().externalPotential().atomEnergy(*atomPtr);
         }
      }
      #endif

      // Reverse sequence of atom types
      for (i = 0; i < nAtom; ++i) {
         assert(molPtr->atom(i).typeId() == atomTypeIds_[i]);
         molPtr->atom(i).setTypeId(atomTypeIds_[nAtom - 1 - i]);   
      }

      // Calculate new energy (with reversed sequence of atom types).
      newEnergy = 0.0;
      #ifndef SIMP_NOPAIR 
      newEnergy += system().pairPotential().moleculeEnergy(*molPtr);
      #endif
      #ifdef SIMP_EXTERNAL
      if (system().hasExternalPotential()) {
         for (i = 0; i < nAtom; ++i) {
            molPtr->atom(i).setTypeId(atomTypeIds_[nAtom - 1 - i]);   
            atomPtr = &molPtr->atom(i);
            newEnergy += system().externalPotential().atomEnergy(*atomPtr);
         }
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
            #ifndef SIMP_NOPAIR
            system().pairPotential().updateAtomCell(*atomPtr);
            #endif
         }

         incrementNAccept();

      }

      return accept;
   }

}
