/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "CfbDoubleRebridgeMove.h"
#include <mcMd/mcSimulation/McSystem.h>
#include <mcMd/simulation/Simulation.h>
#ifndef SIMP_NOPAIR
#include <mcMd/potentials/pair/McPairPotential.h>
#endif
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Bond.h>
#include <mcMd/chemistry/Atom.h>

#include <simp/species/Linear.h>

#include <util/boundary/Boundary.h>
#include <util/global.h>

namespace McMd
{

   using namespace Util;
   using namespace Simp;

   /* 
   * Constructor
   */
   CfbDoubleRebridgeMove::CfbDoubleRebridgeMove(McSystem& system) : 
      CfbRebridgeBase(system),
      speciesId_(-1),
      nRegrow_(-1),
      bridgeLength_(0)
   {  setClassName("CfbDoubleRebridgeMove"); } 
   
   /* 
   * Read parameters speciesId, nRegrow, and nTrial
   */
   void CfbDoubleRebridgeMove::readParameters(std::istream& in) 
   {
      // Read parameters
      readProbability(in);
      read<int>(in, "speciesId", speciesId_);
      read<int>(in, "nRegrow", nRegrow_);
      read<double>(in, "bridgeLength", bridgeLength_);

      // Read parameters hold by RebridgeBase
      CfbRebridgeBase::readParameters(in);

      // Validate
      Linear* chainPtr;
      chainPtr = dynamic_cast<Linear*>(&(simulation().species(speciesId_)));
      if (!chainPtr) {
         UTIL_THROW("Not a Linear species");
      }
  
      // Initialize tables for spring constants and normalizations.
      setup();

      // Allocate array to store old positions 
      iOldPos_.allocate(nRegrow_); 
      jOldPos_.allocate(nRegrow_); 
   }

   /* 
   * Load state from archive.
   */
   void CfbDoubleRebridgeMove::loadParameters(Serializable::IArchive& ar) 
   {
      // Read parameters
      McMove::loadParameters(ar);
      loadParameter<int>(ar, "speciesId", speciesId_);
      loadParameter<int>(ar, "nRegrow", nRegrow_);
      loadParameter<double>(ar, "bridgeLength", bridgeLength_);
      CfbRebridgeBase::loadParameters(ar);

      // Validate
      Linear* chainPtr;
      chainPtr = dynamic_cast<Linear*>(&(simulation().species(speciesId_)));
      if (!chainPtr) {
         UTIL_THROW("Species is not a subclass of Linear");
      }
  
      // Initialize tables for spring constants and normalizations.
      setup();

      // Allocate arrays to store old positions 
      iOldPos_.allocate(nRegrow_); 
      jOldPos_.allocate(nRegrow_); 
   }

   /* 
   * Save state to archive.
   */
   void CfbDoubleRebridgeMove::save(Serializable::OArchive& ar) 
   {
      McMove::save(ar);
      ar & speciesId_;
      ar & nRegrow_;
      ar & bridgeLength_;
      CfbRebridgeBase::save(ar);
   }


   /* 
   * Generate, attempt and accept or reject a Monte Carlo move.
   *
   * A complete move involve the following stepts:
   *
   *    1) choose a molecule-i randomly;
   *    2) find molecule-j of the same type as i, and monomer pairs satisfying
   *       the distance criterion exist;
   *    3) Delete the two rebridging monomers; first molecule-i, then -j.
   *    4) Rebuild molecule-j, then -i. (in the JCP paper, the order is chosen
   *       at random.)
   *
   */
   bool CfbDoubleRebridgeMove::move() 
   {
      bool      found, accept;
      int       iMol, jMol, sign, beginId, nEnd, i;
      Molecule *iPtr;     // pointer to i-th molecule
      Molecule *jPtr;     // pointer to j-th molecule
      int      *bonds;    // bond types of a consequtive blocks of atoms
      Atom     *thisPtr;  // pointer to the current end atom
      Atom     *tempPtr;  // pointer used for swap atom positions
      double    rosenbluth, rosen_r,  rosen_f;
      double    energy, energy_r, energy_f;
      double    po2n, pn2o;
      Vector    swapPos;

      incrementNAttempt();

      // Choose a random direction
      if (random().uniform(0.0, 1.0) > 0.5)
         sign = +1;
      else
         sign = -1;

      // Search for a pair of candidate sites and calculate probability 
      found = forwardScan(sign, iMol, jMol, beginId, po2n);
      if (!found) return found;

      iPtr = &(system().molecule(speciesId_, iMol));
      jPtr = &(system().molecule(speciesId_, jMol));

      // Save positions for atoms to be regrown
      thisPtr = &(iPtr->atom(beginId));
      tempPtr = &(jPtr->atom(beginId));
      for (i = 0; i < nRegrow_; ++i) {
         iOldPos_[i] = thisPtr->position();
         jOldPos_[i] = tempPtr->position();

         thisPtr += sign;
         tempPtr += sign;
      }

      // Prepare for the rebridge move
      accept = false;

      bonds = new int[nRegrow_ + 1];
      if (sign == +1) {
         for (i = 0; i < nRegrow_ + 1; ++i)
            bonds[i] = iPtr->bond(beginId - 1 + nRegrow_ - i).typeId();
      } else {
         for (i = 0; i < nRegrow_ + 1; ++i)
            bonds[i] = iPtr->bond(beginId - nRegrow_ + i).typeId();
      }

      // Deleting atoms: first i-th molecule, then j-th
      rosen_r = 1.0;
      energy_r = 0.0;

      thisPtr = &(iPtr->atom(beginId + (nRegrow_-1)*sign));
      deleteSequence(nRegrow_, sign, thisPtr, bonds, rosenbluth, energy);
      rosen_r *= rosenbluth;
      energy_r += energy;

      thisPtr = &(jPtr->atom(beginId + (nRegrow_-1)*sign));
      deleteSequence(nRegrow_, sign, thisPtr, bonds, rosenbluth, energy);
      rosen_r *= rosenbluth;
      energy_r += energy;

      // Swap the atom positions of the dangling end
      if (sign == +1)
         nEnd = iPtr->nAtom() - beginId - nRegrow_;
      else
         nEnd = beginId + 1 - nRegrow_;

      thisPtr = &(iPtr->atom(beginId + nRegrow_*sign));
      tempPtr = &(jPtr->atom(beginId + nRegrow_*sign));
      for (i = 0; i < nEnd; ++i) {
         swapPos = thisPtr->position();
         //system().moveAtom(*thisPtr, tempPtr->position());
         thisPtr->position() = tempPtr->position();
         #ifndef SIMP_NOPAIR
         system().pairPotential().updateAtomCell(*thisPtr);
         #endif
         //system().moveAtom(*tempPtr, swapPos);
         tempPtr->position() = swapPos;
         #ifndef SIMP_NOPAIR
         system().pairPotential().updateAtomCell(*tempPtr);
         #endif

         thisPtr += sign;
         tempPtr += sign;
      }

      // Regrow atoms: first j-th molecule, then i-th
      rosen_f = 1.0;
      energy_f = 0.0;

      thisPtr = &(jPtr->atom(beginId));
      addSequence(nRegrow_, sign, thisPtr, bonds, rosenbluth, energy);
      rosen_f *= rosenbluth;
      energy_f += energy;

      thisPtr = &(iPtr->atom(beginId));
      addSequence(nRegrow_, sign, thisPtr, bonds, rosenbluth, energy);
      rosen_f *= rosenbluth;
      energy_f += energy;

      // Reverse scan
      found = reverseScan(sign, iMol, jMol, beginId - sign, pn2o);

      // Decide whether to accept or reject
      if (found)
         accept = random().metropolis(rosen_f * pn2o / rosen_r / po2n);
      else
         accept = false;

      if (accept) {

         // Increment counter for accepted moves of this class.
         incrementNAccept();

      } else {

         // If the move is rejected, restore nRegrow_ positions
         thisPtr = &(iPtr->atom(beginId));
         tempPtr = &(jPtr->atom(beginId));
         for (i = 0; i < nRegrow_; ++i) {
            //system().moveAtom(*thisPtr, iOldPos_[i]);
            thisPtr->position() = iOldPos_[i];
            #ifndef SIMP_NOPAIR
            system().pairPotential().updateAtomCell(*thisPtr);
            #endif
            //system().moveAtom(*tempPtr, jOldPos_[i]);
            tempPtr->position() = jOldPos_[i];
            #ifndef SIMP_NOPAIR
            system().pairPotential().updateAtomCell(*tempPtr);
            #endif

            thisPtr += sign;
            tempPtr += sign;
         }

         // Restore atom positions at the dangling end
         thisPtr = &(iPtr->atom(beginId + nRegrow_*sign));
         tempPtr = &(jPtr->atom(beginId + nRegrow_*sign));
         for (i = 0; i < nEnd; ++i) {
            swapPos = thisPtr->position();
            //system().moveAtom(*thisPtr, tempPtr->position());
            thisPtr->position() = tempPtr->position();
            #ifndef SIMP_NOPAIR
            system().pairPotential().updateAtomCell(*thisPtr);
            #endif
            //system().moveAtom(*tempPtr, swapPos);
            tempPtr->position() = swapPos;
            #ifndef SIMP_NOPAIR
            system().pairPotential().updateAtomCell(*tempPtr);
            #endif

            thisPtr += sign;
            tempPtr += sign;
         }

      }

      // Release memory
      delete [] bonds;

      return accept;
   }

   /*
   * Scanning chains to find all possible pairs of rebriding sites satisfying
   * the distance criterion, then choose one randomly to commit the double
   * rebridging move.
   *
   *    probability = 1 / (number of pairs satisfying the distance criterion)
   *
   *    rebriding criterion:
   *
   *          distance(i1,j2) < bridgeLength
   *       && distance(i2,j1) < bridgeLength
   *
   *                   i1   i2
   *           -->------o-o-o----   i-chain
   *                     \ /
   *                      X
   *                     / \
   *           -->------o-o-o----   j-chain
   *                   j1   j2
   *
   *
   */
   bool CfbDoubleRebridgeMove::forwardScan(int sign, int &iMol,
                            int &jMol, int &beginId, double &prob)
   {
      Molecule *iPtr, *jPtr;
      Vector   iPos1, jPos1, iPos2, jPos2;
      int      nMol, nAtom, maxPair;
      int      signRegrow, nPairs, k, choice;
      int      *molBuf, *atomBuf;
      double   lengthSq, bridgeLengthSq;  
      bool     found;

      // Number of molecules
      nMol = system().nMolecule(speciesId_);

      // Randomly choose one molecule
      iPtr = &(system().randomMolecule(speciesId_));
      iMol = system().moleculeId(*iPtr);

      // Number of atoms per molecule
      nAtom = iPtr->nAtom();

      // Require that chain length - 2 >= nRegrow_
      if (nRegrow_ > nAtom - 2) {
         UTIL_THROW("nRegrow_  >= chain length");
      }

      // Allocate arrays holding candidate atom pairs
      maxPair = (nMol - 1) * (nAtom - 1 - nRegrow_);
      molBuf  = new int[maxPair];  
      atomBuf = new int[maxPair];  
      
      // Identify the direction of scaning
      beginId = (sign == -1) ? (nAtom-1) : 0;

      // Auxiliary variables
      signRegrow = sign * (nRegrow_ + 1);
      bridgeLengthSq = bridgeLength_ * bridgeLength_;

      // Nullify counters
      nPairs = 0; 
      found  = false;

      // Loop over j molecules to scan possible pairs
      for (jMol = 0; jMol < nMol; ++jMol) {

         jPtr = &(system().molecule(speciesId_, jMol));

         if (jMol != iMol) {

            for (k = 0; k <  nAtom - 1 - nRegrow_; ++k) {  

               // Check if i1 to j2 bridging is possible
               iPos1 = iPtr->atom(beginId + k*sign).position();
               jPos2 = jPtr->atom(beginId + signRegrow + k*sign).position();
               lengthSq = boundary().distanceSq(iPos1, jPos2);

               if (lengthSq < bridgeLengthSq) {
                  // Check if j2 to i1 bridging is possible          
                  jPos1 = jPtr->atom(beginId + k*sign).position();
                  iPos2 = iPtr->atom(beginId + signRegrow + k*sign).position();
                  lengthSq = boundary().distanceSq(iPos2, jPos1);
                  
                  if (lengthSq < bridgeLengthSq) {
                     molBuf[nPairs]  = jMol;
                     atomBuf[nPairs] = beginId + k*sign;
                     nPairs += 1;
                  }  
               }

            } // loop over atom positions
         } // if jMol != iMol
      } // loop over j molecules
      
      if (nPairs > 0) found = true;
      
      if (found) {
         // iMol is not changed after the random selection
         choice  = random().uniformInt(0, nPairs);  
         jMol    = molBuf[choice];
         beginId = atomBuf[choice] + sign; // return the the active atom ID
         prob    = 1.0 / double(nPairs);
      } else {
         iMol    = -1;  // invalid value
         jMol    = -1;  // invalid value
         beginId = -1;  // invalid value
         prob    = 0.0;
      }   
      
      // Deallocate memory and quit
      delete [] molBuf;
      delete [] atomBuf;
      return found;

   }
   
   /*
   * Scanning chains to find the probability of choosing the specific reverse
   * move.
   *    probability = 1 / (number of pairs satisfying the distance criterion)
   *
   * This method is different than the forward one in two aspects:
   *    (1) It first check if the old configuration is one possible reverse
   *        move; and return "0" probability if not.
   *    (2) It does not store the indices of possible monomer pairs and simply
   *        count the number of all possible bridging sites.
   */
   bool CfbDoubleRebridgeMove::reverseScan(int sign, int iMol,
                                int jMol, int beginId, double &prob)
   {
      Molecule *iPtr, *jPtr;
      Vector   iPos1, jPos1, iPos2, jPos2;
      int      nMol, nAtom, j;
      int      signRegrow, nPairs, headId, k;
      double   lengthSq, bridgeLengthSq;  
      bool     found;

      // Retrieve molecule pointers
      iPtr = &(system().molecule(speciesId_, iMol));
      jPtr = &(system().molecule(speciesId_, jMol));

      // Number of molecules
      nMol = system().nMolecule(speciesId_);

      // Number of atoms per molecule
      nAtom = iPtr->nAtom();

      // Auxiliary variables
      signRegrow = sign * (nRegrow_ + 1);
      bridgeLengthSq = bridgeLength_ * bridgeLength_;

      // Initialize returning variable
      found = false;

      // Check if i1 to j2 bridging is possible
      iPos1 = iPtr->atom(beginId).position();
      jPos2 = jPtr->atom(beginId + signRegrow).position();
      lengthSq = boundary().distanceSq(iPos1, jPos2);

      if (lengthSq < bridgeLengthSq) {
         // Check if j2 to i1 bridging is possible          
         jPos1 = jPtr->atom(beginId).position();
         iPos2 = iPtr->atom(beginId + signRegrow).position();
         lengthSq = boundary().distanceSq(iPos2, jPos1);
                  
         if (lengthSq < bridgeLengthSq) {
            found = true;
         } else {
            prob = 0.0;
            return found;
         }
      } else {
         prob = 0.0;
         return found;
      }

      // Auxiliary variables
      headId = (sign == -1) ? (nAtom-1) : 0;
      iMol = system().moleculeId(*iPtr);

      // Nullify counters
      nPairs = 0; 

      // Loop over j molecules to scan possible pairs
      for (j = 0; j < nMol; ++j) {

         jPtr = &(system().molecule(speciesId_, j));

         if (j != iMol) {

            for (k = 0; k <  nAtom - 1 - nRegrow_; ++k) {  

               // Check if i1 to j2 bridging is possible
               iPos1 = iPtr->atom(headId + k*sign).position();
               jPos2 = jPtr->atom(headId + signRegrow + k*sign).position();
               lengthSq = boundary().distanceSq(iPos1, jPos2);

               if (lengthSq < bridgeLengthSq) {
                  // Check if j2 to j1 bridging is possible          
                  jPos1 = jPtr->atom(headId + k*sign).position();
                  iPos2 = iPtr->atom(headId + signRegrow + k*sign).position();
                  lengthSq = boundary().distanceSq(iPos2, jPos1);
                  
                  if (lengthSq < bridgeLengthSq) nPairs += 1;
               }

            } // loop over atom positions
         } // if j != iMol
      } // loop over j molecules
 
      // nPairs is always greater than or equal to 1
      prob = 1.0 / double(nPairs);
      return found;
   }
   
}
