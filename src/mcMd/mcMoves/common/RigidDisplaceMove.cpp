#ifndef MCMD_RIGID_DISPLACE_MOVE_CPP
#define MCMD_RIGID_DISPLACE_MOVE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "RigidDisplaceMove.h"
#include <mcMd/simulation/Simulation.h>
#include <mcMd/mcSimulation/McSystem.h>
#include <mcMd/mcSimulation/mc_potentials.h>
#include <mcMd/species/Species.h>
#include <util/boundary/Boundary.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <util/space/Dimension.h>
#include <util/global.h>

namespace McMd
{

   using namespace Util;

   /* 
   * Constructor
   */
   RigidDisplaceMove::RigidDisplaceMove(McSystem& system) 
    : SystemMove(system),
      delta_(0.0),
      speciesId_(-1)
   {  setClassName("RigidDisplaceMove"); }
   
   /* 
   * Read speciesId and delta.
   */
   void RigidDisplaceMove::readParameters(std::istream& in) 
   {
      readProbability(in);
      read<int>(in, "speciesId", speciesId_);
      read<double>(in, "delta", delta_);
      readBlank(in);

      nAtom_ = simulation().species(speciesId_).nAtom();
      oldPositions_.allocate(nAtom_);
   }
   
   /* 
   * Generate, attempt and accept or reject a rigid molecule translation.
   */
   bool RigidDisplaceMove::move() 
   { 
      Vector    oldPos, newPos, dr;
      double    newEnergy, oldEnergy;
      Molecule* molPtr;
      Atom*     atomPtr;
      int       iAtom, j;

      incrementNAttempt();

      // Choose a molecule at random
      molPtr = &(system().randomMolecule(speciesId_));

      // Calculate current energy and store old positions.
      // Ignore covalent bonding energy, because it won't change.
      oldEnergy = 0.0;
      for (iAtom = 0; iAtom < nAtom_; ++iAtom) {
         atomPtr = &molPtr->atom(iAtom);
         oldPositions_[iAtom] = atomPtr->position();
         #ifndef INTER_NOPAIR
         oldEnergy += system().pairPotential().atomEnergy(*atomPtr);
         #endif
         #ifdef INTER_EXTERNAL
         oldEnergy += system().externalPotential().atomEnergy(*atomPtr);
         #endif
         #ifdef INTER_TETHER
         oldEnergy += system().atomTetherEnergy(*atomPtr);
         #endif
      }

      // Generate trial displacement Vector dr
      for (j = 0; j < Dimension; ++j) {
         dr[j] = random().uniform(-delta_, delta_);
      }

      // Move every atom by dr and calculate new trial energy.
      newEnergy = 0.0;
      for (iAtom = 0; iAtom < nAtom_; ++iAtom) {
         atomPtr = &molPtr->atom(iAtom);
         atomPtr->position() += dr;
         boundary().shift(atomPtr->position());
         #ifndef INTER_NOPAIR
         newEnergy += system().pairPotential().atomEnergy(*atomPtr);
         #endif
         #ifdef INTER_EXTERNAL
         newEnergy += system().externalPotential().atomEnergy(*atomPtr);
         #endif
         #ifdef INTER_TETHER
         newEnergy += system().atomTetherEnergy(*atomPtr);
         #endif
      }

      // Decide whether to accept the move
      bool accept = random().metropolis(boltzmann(newEnergy - oldEnergy));

      if (accept) {
   
         #ifndef INTER_NOPAIR
         // Update cell
         for (iAtom = 0; iAtom < nAtom_; ++iAtom) {
            system().pairPotential().updateAtomCell(molPtr->atom(iAtom));
         }
         #endif

         incrementNAccept();

      } else {
   
         // Return atoms to original positions
         for (iAtom = 0; iAtom < nAtom_; ++iAtom) {
            molPtr->atom(iAtom).position() = oldPositions_[iAtom];
         }
   
      }
      return accept;
   }

}      
#endif
