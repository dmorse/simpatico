/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "CfbHomoReptationMove.h"
#include <mcMd/mcSimulation/McSystem.h>
#include <mcMd/simulation/Simulation.h>
#ifndef SIMP_NOPAIR
#include <mcMd/potentials/pair/McPairPotential.h>
#endif
#include <mcMd/species/Homopolymer.h>
#include <util/boundary/Boundary.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <mcMd/chemistry/Bond.h>

#include <util/space/Vector.h>
#include <util/global.h>

namespace McMd
{

   using namespace Util;

   /* 
   * Constructor
   */
   CfbHomoReptationMove::CfbHomoReptationMove(McSystem& system) : 
      CfbEndBase(system),
      speciesId_(-1)
   {  setClassName("CfbHomoReptationMove"); } 
   
   /* 
   * Read parameters speciesId, nRegrow, and nTrial
   */
   void CfbHomoReptationMove::readParameters(std::istream& in) 
   {
      // Read parameters
      readProbability(in);
      read<int>(in, "speciesId", speciesId_);
      read<int>(in, "nTrial", nTrial_);

      // Check that 0 < nTrial_ <= MaxTrial
      if (nTrial_ <=0 || nTrial_ > MaxTrial_) {
         UTIL_THROW("Invalid value input for nTrial");
      }

      // Check that the Species is actually a Homopolymer
      Homopolymer* speciesPtr;
      speciesPtr = dynamic_cast<Homopolymer*>(&(simulation().species(speciesId_)));
      if (!speciesPtr) {
         UTIL_THROW("Attempted use of CfbHomoReptation with Species that is not a Homopolymer");
      }
  
   }

   /* 
   * Generate, attempt and accept or reject a Monte Carlo move.
   */
   bool CfbHomoReptationMove::move() 
   {
      Vector    oldPos, newPos;
      double    rosen_r,  rosen_f;
      double    energy_r, energy_f;
      Atom     *tailPtr; // pointer to the tail atom (to be removed)
      Atom     *atomPtr; // resetable atom pointer
      Molecule *molPtr;  // pointer to randomly chosen molecule
      int       length, sign, headId, tailId, bondType, i;
      bool      accept;
     
      incrementNAttempt();

      // Choose a molecule at random
      molPtr = &(system().randomMolecule(speciesId_));
      length = molPtr->nAtom();

      // Choose which chain end to regrow
      if (random().uniform(0.0, 1.0) > 0.5) {
         sign = +1;
         headId = length - 1;
         tailId = 0;
      } else {
         sign = -1;
         headId = 0;
         tailId = length - 1;
      }

      // Store current positions of tail
      oldPos = molPtr->atom(tailId).position();
   
      // Delete tail monomers
      tailPtr = &(molPtr->atom(tailId));
      atomPtr = tailPtr + sign;
      if (sign == 1) {
         bondType = molPtr->bond(0).typeId();
      } else {
         bondType = molPtr->bond(length-2).typeId();
      }
      deleteEndAtom(tailPtr, atomPtr, bondType, rosen_r, energy_r);

      #ifndef SIMP_NOPAIR
      // Delete from McSystem cell list
      system().pairPotential().deleteAtom(*tailPtr);
      #endif
   
      // Regrow head, using tailPtr to store properties of the new head.
      atomPtr = &(molPtr->atom(headId));
      tailPtr->setTypeId(atomPtr->typeId());
      tailPtr->mask().clear();
      #ifndef MCMD_NOMASKBONDED
      tailPtr->mask().append(*atomPtr);
      #endif
      if (sign == 1) {
         bondType = molPtr->bond(length-2).typeId();
      } else {
         bondType = molPtr->bond(0).typeId();
      }
      addEndAtom(tailPtr, atomPtr, bondType, rosen_f, energy_f);

      // Restore original type and connectivity of tail Atom
      atomPtr = tailPtr + sign;
      tailPtr->setTypeId(atomPtr->typeId());
      tailPtr->mask().clear();
      #ifndef MCMD_NOMASKBONDED
      tailPtr->mask().append(*atomPtr);
      #endif

      // Decide whether to accept or reject
      accept = random().metropolis(rosen_f/rosen_r);
      if (accept) {

         // Increment nAccept, the number accepted moves.
         incrementNAccept();

         // Store new head position
         newPos = tailPtr->position();
   
         tailPtr->position() = atomPtr->position();

         #ifndef SIMP_NOPAIR
         // Add back to system cell list
         system().pairPotential().addAtom(*tailPtr);
         #endif
   
         // Shift atom positions towards the head
         for (i=1; i < length - 1; ++i) {
            //system().moveAtom(*atomPtr, (atomPtr+sign)->position());
            atomPtr->position() = (atomPtr+sign)->position();
            #ifndef SIMP_NOPAIR
            system().pairPotential().updateAtomCell(*atomPtr);
            #endif
            atomPtr += sign;
         }
   
         // Move head atom to new chosen position
         //system().moveAtom(*atomPtr, newPos);
         atomPtr->position() = newPos;
         #ifndef SIMP_NOPAIR
         system().pairPotential().updateAtomCell(*atomPtr);
         #endif

      } else {

         // Restore old position of tail
         tailPtr->position() = oldPos;

         #ifndef SIMP_NOPAIR
         // Add tail back to System cell list.
         system().pairPotential().addAtom(*tailPtr);
         #endif
   
      }

      return accept;
   
   }
   
}
