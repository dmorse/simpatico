#ifndef MCMD_CFB_REBRIDGE_MOVE_CPP
#define MCMD_CFB_REBRIDGE_MOVE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, Jian Qin and David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "CfbRebridgeMove.h"
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
   CfbRebridgeMove::CfbRebridgeMove(McSystem& system) : 
      CfbRebridgeBase(system),
      speciesId_(-1),
      nRegrow_(-1)
   {} 
   
   /* 
   * Read parameters speciesId, nRegrow, and nTrial
   */
   void CfbRebridgeMove::readParam(std::istream& in) 
   {
      // Read parameters
      readProbability(in);
      read<int>(in, "speciesId", speciesId_);
      read<int>(in, "nRegrow", nRegrow_);

      // Read parameters hold by RebridgeBase
      CfbRebridgeBase::readParam(in);

      // Initialize tables for spring constants and normalizations.
      setup();

      // Use dynamic_cast to check that the Species is actually a Linear species
      Linear* chainPtr;
      chainPtr = dynamic_cast<Linear*>(&(simulation().species(speciesId_)));
  
      // Allocate array to store old positions 
      oldPos_.allocate(nRegrow_); 
   }


   /* 
   * Generate, attempt and accept or reject a Monte Carlo move.
   *
   * Convention for index
   *       head of molecule                tail of molecule
   *          |                                   |
   *          0   1   2   3 ...                length-1
   *          o---o---o---o---o---o---o---o---o---o
   *              |_______________|___________|
   *                valid beginId    nRegrow
   *
   */
   bool CfbRebridgeMove::move() 
   {
      Molecule *molPtr;   // pointer to randomly chosen molecule
      Atom     *thisPtr;  // pointer to the current end atom
      int       length, sign, beginId, endId, i;
      int       *bonds;
      double    rosen_r,  rosen_f;
      double    energy_r, energy_f;
      bool      accept;

      incrementNAttempt();

      // Choose a molecule at random
      molPtr = &(system().randomMolecule(speciesId_));
      length = molPtr->nAtom();

      // Require that chain length - 2 >= nRegrow_
      if (nRegrow_ > length - 2) {
         UTIL_THROW("nRegrow_  >= chain length");
      }

      // Allocate bonds array
      bonds = new int[nRegrow_ + 1];

      // Choose the beginId of the growing interior bridge
      beginId = random().uniformInt(1, length - nRegrow_);

      // Choose direction: sign = +1 if beginId < endId
      if (random().uniform(0.0, 1.0) > 0.5) {
         sign = +1;
         endId = beginId + (nRegrow_ - 1);
      } else {
         sign = -1;
         endId   = beginId;
         beginId = endId + (nRegrow_ - 1);
      }

      // Store current atomic positions from segment to be regrown
      thisPtr = &(molPtr->atom(beginId));
      for (i = 0; i < nRegrow_; ++i) {
         oldPos_[i] = thisPtr->position();
         thisPtr += sign;
      }

      // Get bond types.
      if (sign == +1) {
         for (i = 0; i < nRegrow_ + 1; ++i)
            bonds[i] = molPtr->bond(endId - i).typeId();
      } else {
         for (i = 0; i < nRegrow_ + 1; ++i)
            bonds[i] = molPtr->bond(endId - 1 + i).typeId();
      }

      // Delete atoms: endId -> beginId
      thisPtr = &(molPtr->atom(endId));
      deleteSequence(nRegrow_, sign, thisPtr, bonds, rosen_r, energy_r);

      // Regrow atoms: endId -> beginId
      thisPtr = &(molPtr->atom(beginId));
      addSequence(nRegrow_, sign, thisPtr, bonds, rosen_f, energy_f);

      // Release bonds array
      delete [] bonds;

      // Decide whether to accept or reject
      accept = random().metropolis(rosen_f/rosen_r);
      if (accept) {

         // Increment counter for accepted moves of this class.
         incrementNAccept();

         // If the move is accepted, keep current positions.

      } else {

         // If the move is rejected, restore old positions
         thisPtr = &(molPtr->atom(beginId));
         for (i = 0; i < nRegrow_; ++i) {
            //system().moveAtom(*thisPtr, oldPos_[i]);
            thisPtr->position() = oldPos_[i];
            #ifndef INTER_NOPAIR
            system().pairPotential().updateAtomCell(*thisPtr);
            #endif
            thisPtr += sign;
         }
   
      }

      return accept;
   }
}
#endif
