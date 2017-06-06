/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "CfbRingRebridgeMove.h"
#include <mcMd/mcSimulation/McSystem.h>
#include <mcMd/simulation/Simulation.h>
#ifndef SIMP_NOPAIR
#include <mcMd/potentials/pair/McPairPotential.h>
#endif
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Bond.h>
#include <mcMd/chemistry/Atom.h>
#include <simp/species/Ring.h>
#include <util/boundary/Boundary.h>
#include <util/global.h>

namespace McMd
{

   using namespace Util;
   using namespace Simp;

   /* 
   * Constructor
   */
   CfbRingRebridgeMove::CfbRingRebridgeMove(McSystem& system) : 
      CfbRebridgeBase(system),
      speciesId_(-1),
      nRegrow_(-1)
   {  setClassName("CfbRingRebridgeMove"); } 
   
   /* 
   * Read parameters speciesId, nRegrow, and nTrial
   */
   void CfbRingRebridgeMove::readParameters(std::istream& in) 
   {
      // Read parameters
      readProbability(in);
      read<int>(in, "speciesId", speciesId_);
      read<int>(in, "nRegrow", nRegrow_);

      // Read parameters hold by RebridgeBase
      CfbRebridgeBase::readParameters(in);

      // Initialize tables for spring constants and normalizations.
      setup();

      // Use dynamic_cast to check that the Species is actually a Linear species
      Ring* ringPtr;
      ringPtr = dynamic_cast<Ring*>(&(simulation().species(speciesId_)));
      if (!ringPtr) {
         UTIL_THROW("Not a Ring species");
      }
  
      // Allocate array to store old positions 
      oldPos_.allocate(nRegrow_); 
   }

   /*
   * Load state from an archive.
   */
   void CfbRingRebridgeMove::loadParameters(Serializable::IArchive& ar)
   {  
      McMove::loadParameters(ar);
      loadParameter<int>(ar, "speciesId", speciesId_);
      loadParameter<int>(ar, "nRegrow", nRegrow_);
      CfbRebridgeBase::loadParameters(ar);

      // Validate
      Ring* ringPtr;
      ringPtr = dynamic_cast<Ring*>(&(simulation().species(speciesId_)));
      if (!ringPtr) {
         UTIL_THROW("Species is not a Ring species");
      }
  
      // Initialize tables for spring constants and normalizations.
      setup();
      oldPos_.allocate(nRegrow_); 
   }

   /*
   * Save state to an archive.
   */
   void CfbRingRebridgeMove::save(Serializable::OArchive& ar)
   {  
      McMove::save(ar);
      ar & speciesId_;
      ar & nRegrow_;
      CfbRebridgeBase::save(ar);
   }

   /* 
   * Generate, attempt and accept or reject a Monte Carlo move.
   *
   * Convention for index
   *                          nAtom
   *          _____________________________________
   *          |                                   |
   *          0   1   2   3 ...                nAtom-1
   *          o---o---o---o---o---o---o---o---o---o
   *  bonds:    0   1   2   3                nAtom-1
   *
   * The algorithm is slightly different from that for linear molecules since a
   * ring molecule does not have ends. The operation of moving the atom pointers
   * have to account for the periodic boundary conditions, which is achived with
   * the help of modId(int s, int n) method.
   *
   */
   bool CfbRingRebridgeMove::move() 
   {
      Molecule *molPtr;    // pointer to randomly chosen molecule
      Atom     *thisPtr;   // pointer to the current end atom
      Atom     *tempPtr;   // pointer to the current end atom
      int      *bonds, iBond;
      int       nAtom, sign, beginId, endId, i;
      double    rosen_r,  rosen_f;
      double    energy_r, energy_f;
      bool      accept;

      incrementNAttempt();

      // Choose a molecule at random
      molPtr = &(system().randomMolecule(speciesId_));
      nAtom  = molPtr->nAtom();

      // Require that chain length - 2 >= nRegrow_
      if (nRegrow_ > nAtom - 2) {
         UTIL_THROW("nRegrow_  >= chain length");
      }

      // Allocate bonds array
      bonds = new int[nRegrow_ + 1];

      // Choose the beginId of the growing interior bridge
      beginId = random().uniformInt(0, nAtom);

      // Choose direction and fill bonds array: sign = +1 if beginId < endId
      if (random().uniform(0.0, 1.0) > 0.5) {
         sign = +1;
         endId = (beginId + (nRegrow_ - 1)) % nAtom;

         // Save bond types
         iBond = endId;
         for (i = 0; i < nRegrow_ + 1; ++i) {
            bonds[i] = molPtr->bond(iBond).typeId();
            iBond--;
            if (iBond < 0) iBond += nAtom;
         }
      } else {
         sign = -1;
         endId = beginId;
         beginId = (endId + (nRegrow_ - 1)) % nAtom;

         // Save bond types
         iBond = endId - 1;
         if (iBond < 0) iBond += nAtom;
         for (i = 0; i < nRegrow_ + 1; ++i) {
            bonds[i] = molPtr->bond(iBond).typeId();
            iBond++;
            if (iBond >= nAtom) iBond -= nAtom;
         }
      }

      // Store current atomic positions from segment to be regrown
      tempPtr = &(molPtr->atom(beginId));
      for (i = 0; i < nRegrow_; ++i) {
         thisPtr = tempPtr + modId(beginId + i*sign, nAtom) - beginId;
         oldPos_[i] = thisPtr->position();
      }

      // Delete atoms: endId -> beginId
      deleteSequence(sign, molPtr, endId, bonds, rosen_r, energy_r);

      // Regrow atoms: endId -> beginId
      addSequence(sign, molPtr, beginId, bonds, rosen_f, energy_f);

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
         tempPtr = &(molPtr->atom(beginId));
         for (i = 0; i < nRegrow_; ++i) {
            thisPtr = tempPtr + modId(beginId + i*sign, nAtom) - beginId;
            //system().moveAtom(*thisPtr, oldPos_[i]);
            thisPtr->position() =  oldPos_[i];
            #ifndef SIMP_NOPAIR
            system().pairPotential().updateAtomCell(*thisPtr);
            #endif
         }
   
      }

      return accept;
   }

   /* 
   * Delete a consecutive sequence of atoms in a Ring.
   */
   void CfbRingRebridgeMove::deleteSequence(int sign, Molecule *molPtr,
          int endId, int *bonds, double &rosenbluth, double &energy)
   {
      Atom   *endPtr;
      Atom   *thisPtr;
      Atom   *prevPtr;
      Atom   *nextPtr;
      int     prevBType, nextBType, nAtom;
      double  rosen_r, energy_r;

      // Initialize weight
      rosenbluth = 1.0;
      energy = 0.0;

      nAtom  = molPtr->nAtom();
      endPtr = &(molPtr->atom(endId));

      // Step 1: delete the last atom
      if (nRegrow_ >= 1) {
         thisPtr = endPtr;
         prevPtr = endPtr + modId(endId - sign,nAtom) - endId;
         nextPtr = endPtr + modId(endId + sign,nAtom) - endId;

         nextBType = bonds[0];
         prevBType = bonds[1];
 
         // Orientation biased trimer rebridge move
         deleteMiddleAtom(thisPtr, prevPtr, nextPtr,
                prevBType, nextBType, rosen_r, energy_r);
         #ifndef SIMP_NOPAIR
         system().pairPotential().deleteAtom(*thisPtr);
         #endif
         rosenbluth *= rosen_r;
         energy += energy_r;
      }

      // step 2: delete the remaining atoms
      for (int i = 0; i < nRegrow_ - 1; ++i) {
         thisPtr = endPtr + modId(endId - (i+1)*sign, nAtom) - endId;
         prevPtr = endPtr + modId(endId - (i+2)*sign, nAtom) - endId;

         prevBType = bonds[i+2];
         deleteEndAtom(thisPtr, prevPtr, prevBType, rosen_r, energy_r);
         #ifndef SIMP_NOPAIR
         system().pairPotential().deleteAtom(*thisPtr);
         #endif

         rosenbluth *= rosen_r;
         energy += energy_r;
      }
   }

   /* 
   * Add a consecutive sequence of atoms in a Ring.
   */
   void CfbRingRebridgeMove::addSequence(int sign, Molecule *molPtr,
          int beginId, int *bonds, double &rosenbluth, double &energy)
   {
      Atom   *beginPtr;
      Atom   *thisPtr;
      Atom   *prevPtr;
      Atom   *nextPtr;
      int     prevBType, nextBType, nAtom;
      double  rosen_f, energy_f;

      // Initialize weight
      rosenbluth = 1.0;
      energy = 0.0;

      // Step 1: grow the first nRegrow_ - 1 atoms
      nAtom = molPtr->nAtom();
      beginPtr = &(molPtr->atom(beginId));

      for (int i = 0; i < nRegrow_ - 1; ++i) {
         thisPtr = beginPtr + modId(beginId+i*sign,nAtom) - beginId;
         prevPtr = beginPtr + modId(beginId+(i-1)*sign,nAtom) - beginId;

         prevBType = bonds[nRegrow_ - i];
         addEndAtom(thisPtr, prevPtr, prevBType, rosen_f, energy_f);
         #ifndef SIMP_NOPAIR
         system().pairPotential().addAtom(*thisPtr);
         #endif

         rosenbluth *= rosen_f;
         energy += energy_f;
      }

      // Step 2: grow the last atoms
      if (nRegrow_ >= 1) {
         prevPtr = beginPtr + modId(beginId+(nRegrow_-2)*sign,nAtom) - beginId;
         thisPtr = beginPtr + modId(beginId+(nRegrow_-1)*sign,nAtom) - beginId;
         nextPtr = beginPtr + modId(beginId+nRegrow_*sign,nAtom) - beginId;

         prevBType = bonds[1];
         nextBType = bonds[0];

         // Invoke the orientation biased trimer re-bridging move
         addMiddleAtom(thisPtr, prevPtr, nextPtr,
                prevBType, nextBType, rosen_f, energy_f);
         #ifndef SIMP_NOPAIR
         system().pairPotential().addAtom(*thisPtr);
         #endif

         rosenbluth *= rosen_f;
         energy += energy_f;
      }
   }

}
