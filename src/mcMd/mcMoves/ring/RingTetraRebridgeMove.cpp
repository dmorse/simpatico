/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RingTetraRebridgeMove.h"
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

   /* 
   * Constructor
   */
   RingTetraRebridgeMove::RingTetraRebridgeMove(McSystem& system) : 
      GroupRebridgeBase(system),
      speciesId_(-1),
      nAtom_(0),
      upperBridge_(0),
      lowerBridge_(0)
   {  setClassName("RingTetraRebridgeMove"); } 
   
   /* 
   * Read parameters speciesId, and bridgeLength_
   */
   void RingTetraRebridgeMove::readParameters(std::istream& in) 
   {
      // Read parameters
      readProbability(in);
      read<int>(in, "speciesId", speciesId_);
      read<double>(in, "upperBridge", upperBridge_);
      read<double>(in, "lowerBridge", lowerBridge_);

      // Use dynamic_cast to check that the Species is actually a Linear species
      Ring* ringPtr;
      ringPtr = dynamic_cast<Ring*>(&(simulation().species(speciesId_)));
      if (!ringPtr) {
         UTIL_THROW("Not a Ring species");
      }

      // Allocate array for atom position vectors.
      nAtom_ = system().simulation().species(speciesId_).nAtom();
      if (nAtom_ < 4) UTIL_THROW("nAtom < 4 for tetra rebridge move.");
      R_.allocate(nAtom_); 

   }

   /*
   * Load state from an archive.
   */
   void RingTetraRebridgeMove::loadParameters(Serializable::IArchive& ar)
   {  
      McMove::loadParameters(ar);
      loadParameter<int>(ar, "speciesId", speciesId_);
      loadParameter<double>(ar, "upperBridge", upperBridge_);
      loadParameter<double>(ar, "lowerBridge", lowerBridge_);
      ar & nAtom_;

      // Validate
      if (nAtom_ != system().simulation().species(speciesId_).nAtom()) {
         UTIL_THROW("Inconsistent values of nAtom");
      }
      Ring* ringPtr;
      ringPtr = dynamic_cast<Ring*>(&(simulation().species(speciesId_)));
      if (!ringPtr) {
         UTIL_THROW("Species is not a Ring species");
      }
      if (nAtom_ < 4) {
         UTIL_THROW("nAtom < 4 for tetra rebridge move.");
      }

      // Allocate array for atom position vectors.
      R_.allocate(nAtom_); 
   }

   /*
   * Save state to an archive.
   */
   void RingTetraRebridgeMove::save(Serializable::OArchive& ar)
   {
      McMove::save(ar);
      ar & speciesId_;
      ar & upperBridge_;
      ar & lowerBridge_;
      ar & nAtom_;
   }


   /* 
   * Generate, attempt and accept or reject a Monte Carlo move.
   *
   * A complete move involve the following stepts:
   *
   *    1) choose a molecule-i randomly;
   *    2) find molecule-j of the same type as i, and monomer pairs satisfying
   *       the distance criterion exist;
   *
   */
   bool RingTetraRebridgeMove::move() 
   {
      Molecule *molPtr;    // pointer to the rebridging molecule
      Atom     *aPtr, *bPtr, *cPtr, *dPtr; // atom pointers
      Atom     *lPtr, *hPtr, *lTmp, *hTmp;
      Vector   swapV;
      int      ia, ic, ib, id;
      int      ibc, ida, iLow, iHigh, nSwap;
      int      bondType;
      double   energy;
      bool     found, accept;

      incrementNAttempt();

      // Randomly pick up a molecule.
      molPtr = &(system().randomMolecule(speciesId_));

      // Randomly choose the a atom in the tetra-group.
      ia = random().uniformInt(0, nAtom_);

      // Search the rebriding sites.
      found = scanBridge(molPtr, ia, ic);
      if (!found) return false;

      ib = modId(ia+1, nAtom_);
      id = modId(ic+1, nAtom_);

      // Calcuate the energy change.
      aPtr = &(molPtr->atom(ia));
      bPtr = &(molPtr->atom(ib));
      cPtr = &(molPtr->atom(ic));
      dPtr = &(molPtr->atom(id));

      bondType = molPtr->bond(ia).typeId();
      tetraEnergy(aPtr, bPtr, cPtr, dPtr, bondType, energy);

      // Decide whether to accept or reject.
      accept = random().metropolis(boltzmann(energy));

      if (accept) {

         // Find the shorter half loop.
         ibc = ic - ib; 
         if (ibc < 0) ibc += nAtom_;

         ida = ia - id; 
         if (ida < 0) ida += nAtom_;

         if (ibc <= ida) {
            lPtr = bPtr;
            hPtr = cPtr;
            iLow = ib;
            iHigh = ic;
            if (ibc%2 == 0)
               nSwap = ibc / 2;
            else
               nSwap = (ibc + 1) / 2;
         } else {
            lPtr = dPtr;
            hPtr = aPtr;
            iLow = id;
            iHigh = ia;
            if (ida%2 == 0)
               nSwap = ida / 2;
            else
               nSwap = (ida + 1) / 2;
         }

         // Swap positions.
         for (int j = 0; j < nSwap; ++j) {
            // moving pointer close to the atom with higher index 
            hTmp = hPtr + modId(iHigh - j, nAtom_) - iHigh;

            // moving pointer close to the atom with lower index
            lTmp = lPtr + modId(iLow + j, nAtom_) - iLow;

            swapV = lTmp->position();

            //system().moveAtom(*lTmp, hTmp->position());
            lTmp->position() = hTmp->position();
            #ifndef SIMP_NOPAIR
            system().pairPotential().updateAtomCell(*lTmp);
            #endif

            //system().moveAtom(*hTmp, swapV);
            hTmp->position() = swapV;
            #ifndef SIMP_NOPAIR
            system().pairPotential().updateAtomCell(*hTmp);
            #endif

         }

         // Increment counter for accepted moves of this class.
         incrementNAccept();

      }

      return accept;
   }

   /*
   * Scan the valid rebridging sits.
   */
   bool RingTetraRebridgeMove::
   scanBridge(Molecule* molPtr, int ia, int &ic)
   {
      Vector    r1, r2, dR;
      int       ib, id, i;
      int       *idList, nGroup;
      double    lenSq, minSq, maxSq;
      bool      found;

      // Preparation.
      found = false;
      maxSq = upperBridge_ * upperBridge_;
      minSq = lowerBridge_ * lowerBridge_;

      // b Atom id.
      ib = modId(ia+1, nAtom_);
      r1 = molPtr->atom(ia).position();
      r2 = molPtr->atom(ib).position();

      // Check the ab bond.
      lenSq = boundary().distanceSq(r1, r2);
      if (lenSq >= maxSq || lenSq <= minSq) {
         ic = -1;       // invalid value
         return found;
      }

      // Retrace the non-periodic molecule shape.
      R_[0] = molPtr->atom(0).position();
      for (i = 1; i < nAtom_; ++i) {
         r1 = molPtr->atom(i-1).position();
         r2 = molPtr->atom(i).position();
         system().boundary().distanceSq(r2, r1, dR);
         R_[i].add(R_[i-1], dR);
      }

      // Loop over all monomers to find the proper rebridging sites.
      nGroup = 0;
      idList = new int[nAtom_ - 3];
      for (i = 1; i < nAtom_ - 2; ++i) {
         ic = modId(ib+i,nAtom_);
         id = modId(ic+1,nAtom_);

         // Distance test.
         dR.subtract(R_[ia],R_[ic]);
         lenSq = dR.square();
         if (lenSq > minSq && lenSq < maxSq) {

         dR.subtract(R_[id],R_[ic]);
         lenSq = dR.square();
         if (lenSq > minSq && lenSq < maxSq) {

         dR.subtract(R_[id],R_[ib]);
         lenSq = dR.square();
         if (lenSq > minSq && lenSq < maxSq) {

            idList[nGroup] = ic; 
            ++nGroup;

         }
         }
         }

      } // Loop over all possible c atoms.

      if (nGroup > 0) {
         found = true;
         ic = idList[random().uniformInt(0,nGroup)];
      } else {
         found = false;
         ic = -1;         // Invalid value.
      }

      delete [] idList;
      return found;

   }
   
}
