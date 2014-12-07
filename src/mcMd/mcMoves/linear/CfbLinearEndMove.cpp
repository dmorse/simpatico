/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "CfbLinearEndMove.h"
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
   CfbLinearEndMove::CfbLinearEndMove(McSystem& system) : 
      CfbLinear(system),
      nRegrow_(-1)
   {  setClassName("CfbLinearEndMove"); } 
   
   /* 
   * Read parameters speciesId, nRegrow, and nTrial
   */
   void CfbLinearEndMove::readParameters(std::istream& in) 
   {
      // Read parameters
      readProbability(in);
      CfbLinear::readParameters(in);
      read<int>(in, "nRegrow", nRegrow_);

      // Allocate 
      oldPos_.allocate(nRegrow_); 
   }

   /* 
   * Read parameters speciesId, nRegrow, and nTrial
   */
   void CfbLinearEndMove::loadParameters(Serializable::IArchive& ar) 
   {
      // Read parameters
      McMove::loadParameters(ar);
      CfbLinear::loadParameters(ar);
      loadParameter<int>(ar, "nRegrow", nRegrow_);

      // Validate

      // Allocate array to store old positions 
      oldPos_.allocate(nRegrow_); 
   }

   /* 
   * Save state to archive.
   */
   void CfbLinearEndMove::save(Serializable::OArchive& ar) 
   {
      McMove::save(ar);
      CfbLinear::save(ar);
      ar & nRegrow_;
   }

   /* 
   * Generate, attempt and accept or reject a Monte Carlo move.
   */
   bool CfbLinearEndMove::move() 
   {
      double rosenbluth, rosen_r,  rosen_f;
      double energy, energy_r, energy_f;
      Molecule *molPtr; // pointer to randomly chosen molecule
      Atom* atom0Ptr; // pointer to the current end atom
      Atom* atom1Ptr; // pointer to "pivot" atom, which is bonded to the end atom
      int nAtom, sign, beginId, endId, i;
      bool accept;
     
      incrementNAttempt();
    
      // Choose a molecule at random
      molPtr = &(system().randomMolecule(speciesId()));
      nAtom = molPtr->nAtom();

      // Require that chain nAtom > nRegrow
      if (nRegrow_ >= nAtom) {
         UTIL_THROW("nRegrow_  >= chain nAtom");
      }

      // Choose which chain end to regrow
      if (random().uniform(0.0, 1.0) > 0.5) {
         sign = +1;
         beginId = nAtom - nRegrow_;
         endId   = nAtom - 1;
      } else {
         sign = -1;
         beginId = nRegrow_ - 1;
         endId   = 0;
      }
   
      // Store current atomic positions from segment to be regrown
      for (i = beginId; i <= endId; ++i) {
         oldPos_[i] = molPtr->atom(i).position();
      }
   
      // Delete monomers, starting from chain end
      rosen_r  = 1.0;
      energy_r = 0.0;
      for (i = endId; i >= beginId; --i) {
         deleteAtom(*molPtr, i, sign, rosenbluth, energy);
         #ifndef INTER_NOPAIR
         system().pairPotential().deleteAtom(molPtr->atom(i));
         #endif
         rosen_r *= rosenbluth;
         energy_r += energy;
      }
   
      // Regrow monomers
      rosen_f  = 1.0;
      energy_f = 0.0;
      for (i = beginId; i <= endId; ++i) {
         atom0Ptr = &molPtr->atom(i);
         atom1Ptr = atom0Ptr - sign;
         addAtom(*molPtr, *atom0Ptr, *atom1Ptr, i, sign, rosenbluth, energy);
         rosen_f  *= rosenbluth;
         energy_f += energy;

         #ifndef INTER_NOPAIR
         // Add end atom to McSystem cell list
         system().pairPotential().addAtom(*atom0Ptr);
         #endif
      }
   
      // Decide whether to accept or reject
      accept = random().metropolis(rosen_f/rosen_r);
      if (accept) {

         // Increment counter for accepted moves of this class.
         incrementNAccept();

         // If the move is accepted, keep current positions.

      } else {

         // If the move is rejected, restore old positions
         for (i = beginId; i <= endId; ++i) {
            atom0Ptr = &(molPtr->atom(i));
            atom0Ptr->position() = oldPos_[i];
            #ifndef INTER_NOPAIR
            system().pairPotential().updateAtomCell(*atom0Ptr);
            #endif
         }
   
      }

      return accept;
   }
}
