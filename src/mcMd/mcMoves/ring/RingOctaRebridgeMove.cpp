/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RingOctaRebridgeMove.h"
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
   RingOctaRebridgeMove::RingOctaRebridgeMove(McSystem& system) : 
      GroupRebridgeBase(system),
      speciesId_(-1),
      upperBridge_(0),
      lowerBridge_(0)
   {  setClassName("RingOctaRebridgeMove"); } 
   
   /* 
   * Read parameters speciesId, and nTrial
   */
   void RingOctaRebridgeMove::readParameters(std::istream& in) 
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
   }

   /*
   * Load state from an archive.
   */
   void RingOctaRebridgeMove::loadParameters(Serializable::IArchive& ar)
   {  
      McMove::loadParameters(ar);
      loadParameter<int>(ar, "speciesId", speciesId_);
      loadParameter<double>(ar, "upperBridge", upperBridge_);
      loadParameter<double>(ar, "lowerBridge", lowerBridge_);

      // Validate
      Ring* ringPtr;
      ringPtr = dynamic_cast<Ring*>(&(simulation().species(speciesId_)));
      if (!ringPtr) {
         UTIL_THROW("Species is not a Ring species");
      }
   }

   /*
   * Save state to an archive.
   */
   void RingOctaRebridgeMove::save(Serializable::OArchive& ar)
   {
      McMove::save(ar);
      ar & speciesId_;
      ar & upperBridge_;
      ar & lowerBridge_;
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
   bool RingOctaRebridgeMove::move() 
   {
      Molecule *molPtr1;    // pointer to the rebridging molecule
      Molecule *molPtr2;    // pointer to the rebridging molecule
      Atom     *mPtr, *nPtr; // atom pointers
      Vector   swapV;
      int      nAtom, im, in, molId2;
      double   energy;
      bool     found, accept;

      incrementNAttempt();

      // Randomly pick up a molecule.
      molPtr1 = &(system().randomMolecule(speciesId_));
      nAtom = molPtr1->nAtom();

      // Randomly choose the a atom in the tetra-group.
      im = random().uniformInt(0, nAtom);

      // Search the rebriding sites.
      found = scanBridge(molPtr1, im, molId2, in);
      if (!found) return false;

      // Calcuate the energy change.
      mPtr = &(molPtr1->atom(im));
      molPtr2 = &(system().molecule(speciesId_, molId2));
      nPtr = &(molPtr2->atom(in));

      energy = 0.0;
      energy -= system().atomPotentialEnergy(*mPtr);
      energy -= system().atomPotentialEnergy(*nPtr);

      // Swap positions of atoms m and n.
      swapV = mPtr->position();

      //system().moveAtom(*mPtr, nPtr->position());
      mPtr->position() = nPtr->position();
      #ifndef SIMP_NOPAIR
      system().pairPotential().updateAtomCell(*mPtr);
      #endif

      //system().moveAtom(*nPtr, swapV);
      nPtr->position() = swapV;
      #ifndef SIMP_NOPAIR
      system().pairPotential().updateAtomCell(*nPtr);
      #endif

      energy += system().atomPotentialEnergy(*mPtr);
      energy += system().atomPotentialEnergy(*nPtr);

      // Decide whether to accept or reject.
      accept = random().metropolis(boltzmann(energy));

      if (accept) {
         // Increment counter for accepted moves of this class.
         incrementNAccept();
      } else {
         // Swap back bead positions.
         swapV = mPtr->position();

         //system().moveAtom(*mPtr, nPtr->position());
         mPtr->position() = nPtr->position();
         #ifndef SIMP_NOPAIR
         system().pairPotential().updateAtomCell(*mPtr);
         #endif

         //system().moveAtom(*nPtr, swapV);
         nPtr->position() = swapV;
         #ifndef SIMP_NOPAIR
         system().pairPotential().updateAtomCell(*nPtr);
         #endif
      }

      return accept;
   }

   /*
   * Return true if the octahedron formed by the 6 closely approaching atoms
   * satisfy the 8 distance criterion.
   */
   bool RingOctaRebridgeMove::
   scanBridge(Molecule* molPtr1, int im, int& molId2, int &in)
   {
      Molecule *molPtr2;
      Atom     *pn;
      Vector   aPos, bPos, cPos, dPos, mPos, nPos;
      int      nAtom, ia, ib, ic, id, dmn, dnm;
      int      nNeighbor, j, *idList, *molIdList;
      int      molId1, choice, nGroup;
      double   lenSq, minSq, maxSq;
      bool     validIndex, found;

      // Number of atoms per molecule.
      nAtom = molPtr1->nAtom();

      // The a and b Atom id.
      ia = modId(im-1, nAtom);
      ib = modId(im+1, nAtom);

      aPos = molPtr1->atom(ia).position();
      bPos = molPtr1->atom(ib).position();
      mPos = molPtr1->atom(im).position();

      // Prepare for the scan.
      found = false;
      maxSq = upperBridge_ * upperBridge_;
      minSq = lowerBridge_ * lowerBridge_;

      // Check the am bond.
      lenSq = boundary().distanceSq(aPos, mPos);
      if (lenSq >= maxSq || lenSq <= minSq) {
         in = -1;       // invalid value
         molId2 = -1;   // invalid value
         return found;
      }

      // Check the mb bond.
      lenSq = boundary().distanceSq(mPos, bPos);
      if (lenSq >= maxSq || lenSq <= minSq) {
         in = -1;       // invalid value
         molId2 = -1;   // invalid value
         return found;
      }

      // Get the neighbor list.
      #ifndef SIMP_NOPAIR
      system().pairPotential().cellList().getNeighbors(mPos, neighbors_);
      nNeighbor = neighbors_.size();
      idList = new int[nNeighbor];
      molIdList = new int[nNeighbor];

      // Loop over neighboring monomers and identify the octa-group.
      nGroup = 0;
      molId1 = molPtr1->id();
      for (j = 0; j < nNeighbor; ++j) {
         pn = neighbors_[j];
         in = pn->indexInMolecule();
         molId2 = pn->molecule().id();

         validIndex = true;

         if (molId2 == molId1) {
            // Check if im and in are separated by at least two atoms.
            dmn = in - im;
            if (dmn < 0) dmn += nAtom;

            dnm = im - in;
            if (dnm < 0) dnm += nAtom;

            if (dmn <= 2 && dnm <= 2) validIndex = false;
         }

         if (pn->molecule().species().id() != speciesId_)
            validIndex = false;

         if (validIndex) {

            molPtr2 = &(system().molecule(speciesId_, molId2));

            ic = modId(in-1, nAtom);
            id = modId(in+1, nAtom);

            cPos = molPtr2->atom(ic).position();
            dPos = molPtr2->atom(id).position();
            nPos = molPtr2->atom(in).position();

            // Distance test.
            lenSq = boundary().distanceSq(cPos, mPos);
            if (lenSq > minSq && lenSq < maxSq) {

            lenSq = boundary().distanceSq(mPos, dPos);
            if (lenSq > minSq && lenSq < maxSq) {

            lenSq = boundary().distanceSq(aPos, nPos);
            if (lenSq > minSq && lenSq < maxSq) {

            lenSq = boundary().distanceSq(nPos, bPos);
            if (lenSq > minSq && lenSq < maxSq) {

            lenSq = boundary().distanceSq(cPos, nPos);
            if (lenSq > minSq && lenSq < maxSq) {

            lenSq = boundary().distanceSq(nPos, dPos);
            if (lenSq > minSq && lenSq < maxSq) {

               idList[nGroup] = in; 
               molIdList[nGroup] = molId2; 
               ++nGroup;

            } } } } } }  // 6 additional distance test

         } // if the atom id "in" is valid

      } // loop over neighbors
      #endif

      if (nGroup > 0) {
         found = true;
         choice = random().uniformInt(0,nGroup);
         in = idList[choice];
         molId2 = molIdList[choice];
      } else {
         found = false;
         in = -1;         // Invalid value.
         molId2 = -1;     // Invalid value.
      }

      delete [] molIdList;
      delete [] idList;
      return found;

   }
   
}
