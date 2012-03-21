#ifndef CFB_END_MOVE_CPP
#define MCMD_CFB_END_MOVE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "CfbEndMove.h"
#include <mcMd/mcSimulation/McSystem.h>
#include <mcMd/simulation/Simulation.h>
#ifndef INTER_NOPAIR
#include <mcMd/potentials/pair/McPairPotential.h>
#endif
#include <mcMd/species/Linear.h>
#include <util/boundary/Boundary.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Bond.h>
#include <mcMd/chemistry/Atom.h>
#include <util/global.h>

namespace McMd
{

   using namespace Util;

   /* 
   * Constructor
   */
   CfbEndMove::CfbEndMove(McSystem& system) : 
      CfbEndBase(system),
      speciesId_(-1),
      nRegrow_(-1)
   {} 
   
   /* 
   * Read parameters speciesId, nRegrow, and nTrial
   */
   void CfbEndMove::readParam(std::istream& in) 
   {
      // Read parameters
      readProbability(in);
      read<int>(in, "speciesId", speciesId_);
      read<int>(in, "nRegrow", nRegrow_);
      read<int>(in, "nTrial", nTrial_);

      // Check that 0 < nTrial_ <= MaxTrial
      if (nTrial_ <=0 || nTrial_ > MaxTrial_) {
         UTIL_THROW("Invalid value input for nTrial");
      }

      // Use dynamic_cast to check that the Species is actually a Linear species
      Linear* chainPtr;
      chainPtr = dynamic_cast<Linear*>(&(simulation().species(speciesId_)));
  
      // Allocate array to store old positions 
      oldPos_.allocate(nRegrow_); 
   }

   /* 
   * Generate, attempt and accept or reject a Monte Carlo move.
   */
   bool CfbEndMove::move() 
   {
      double    rosenbluth, rosen_r,  rosen_f;
      double    energy, energy_r, energy_f;
      Atom     *endPtr;   // pointer to the current end atom
      Atom     *pvtPtr;   // pointer to "pivot" atom, which is bonded to the end atom
      Molecule *molPtr;   // pointer to randomly chosen molecule
      int       length, sign, beginId, endId, bondType, i;
      bool      accept;
     
      incrementNAttempt();
    
      // Choose a molecule at random
      molPtr = &(system().randomMolecule(speciesId_));
      length = molPtr->nAtom();

      // Require that chain length > nRegrow
      if (nRegrow_ >= length) {
         UTIL_THROW("nRegrow_  >= chain length");
      }

      // Choose which chain end to regrow
      if (random().uniform(0.0, 1.0) > 0.5) {
         sign = +1;
         beginId = length - nRegrow_;
         endId   = length - 1;
      } else {
         sign = -1;
         beginId = nRegrow_ - 1;
         endId   = 0;
      }
   
      // Store current atomic positions from segment to be regrown
      endPtr = &(molPtr->atom(beginId));
      for (i = 0; i < nRegrow_; ++i) {
         oldPos_[i] = endPtr->position();
         endPtr += sign;
      }
   
      // Delete monomers, starting from chain end
      rosen_r  = 1.0;
      energy_r = 0.0;
      endPtr = &(molPtr->atom(endId));
      for (i = 0; i < nRegrow_; ++i) {
         pvtPtr  = endPtr - sign;
         if (sign == 1) {
            bondType = molPtr->bond(endId - 1 - i).typeId();
         } else {
            bondType = molPtr->bond(i).typeId();
         }
         deleteEndAtom(endPtr, pvtPtr, bondType, rosenbluth, energy);
         #ifndef INTER_NOPAIR
         system().pairPotential().deleteAtom(*endPtr);
         #endif
         rosen_r  *= rosenbluth;
         energy_r += energy;
         endPtr -= sign;
      }
   
      // Regrow monomers
      rosen_f  = 1.0;
      energy_f = 0.0;
      endPtr = &(molPtr->atom(beginId));
      for (i = 0; i < nRegrow_; ++i) {
         pvtPtr = endPtr - sign;
         if (sign == 1) {
            bondType = molPtr->bond(endId - 1 - i).typeId();
         } else {
            bondType = molPtr->bond(i).typeId();
         }
         addEndAtom(endPtr, pvtPtr, bondType, rosenbluth, energy);
         rosen_f  *= rosenbluth;
         energy_f += energy;

         #ifndef INTER_NOPAIR
         // Add end atom to McSystem cell list
         system().pairPotential().addAtom(*endPtr);
         #endif

         endPtr += sign;
      }
   
      // Decide whether to accept or reject
      accept = random().metropolis(rosen_f/rosen_r);
      if (accept) {

         // Increment counter for accepted moves of this class.
         incrementNAccept();

         // If the move is accepted, keep current positions.

      } else {

         // If the move is rejected, restore old positions
         endPtr = &(molPtr->atom(beginId));
         for (i = 0; i < nRegrow_; ++i) {
            // system().moveAtom(*endPtr, oldPos_[i]);
            endPtr->position() = oldPos_[i];
            #ifndef INTER_NOPAIR
            system().pairPotential().updateAtomCell(*endPtr);
            #endif
            endPtr += sign;
         }
   
      }

      return accept;
   
   }
}
#endif
