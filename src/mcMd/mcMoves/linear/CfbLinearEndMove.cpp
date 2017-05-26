/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "CfbLinearEndMove.h"
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

      // Validate
      if (nRegrow_ <= 0) {
         UTIL_THROW("nRegrow_  <= 0");
      }
      int nAtom = system().simulation().species(speciesId()).nAtom();
      if (nRegrow_ >= nAtom) {
         UTIL_THROW("nRegrow_  >= nAtom");
      }

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
      if (nRegrow_ <= 0) {
         UTIL_THROW("nRegrow_  <= 0");
      }
      int nAtom = system().simulation().species(speciesId()).nAtom();
      if (nRegrow_ >= nAtom) {
         UTIL_THROW("nRegrow_  >= nAtom");
      }

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
      Atom* atom1Ptr; // pointer to "pivot" atom, which is bonded to atom0
      int nAtom, sign, beginId, endId, atomId, i;
      bool accept;

      incrementNAttempt();

      // Choose a molecule at random
      molPtr = &(system().randomMolecule(speciesId()));
      nAtom = molPtr->nAtom();

      // Require that chain nAtom > nRegrow
      if (nRegrow_ >= nAtom) {
         UTIL_THROW("nRegrow_  >= nAtom");
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
      atomId = beginId;
      for (i = 0; i < nRegrow_; ++i) {
         oldPos_[i] = molPtr->atom(atomId).position();
         atomId += sign;
      }

      // Delete monomers, starting from chain end
      rosen_r = 1.0;
      energy_r = 0.0;
      atomId = endId;
      for (i = 0; i < nRegrow_; ++i) {
         deleteAtom(*molPtr, atomId, sign, rosenbluth, energy);
         #ifndef SIMP_NOPAIR
         // Add end atom from cell list
         system().pairPotential().deleteAtom(molPtr->atom(atomId));
         #endif
         rosen_r *= rosenbluth;
         energy_r += energy;
         atomId -= sign;
      }
      assert(atomId == beginId - sign);

      // Regrow monomers
      rosen_f = 1.0;
      energy_f = 0.0;
      atomId = beginId;
      for (i = 0; i < nRegrow_; ++i) {
         atom0Ptr = &molPtr->atom(atomId);
         atom1Ptr = atom0Ptr - sign;
         addAtom(*molPtr, *atom0Ptr, *atom1Ptr, atomId, sign,
                 rosenbluth, energy);
         rosen_f *= rosenbluth;
         energy_f += energy;
         #ifndef SIMP_NOPAIR
         // Add end atom to cell list
         system().pairPotential().addAtom(*atom0Ptr);
         #endif
         atomId += sign;
      }
      assert(atomId == endId + sign);

      // Accept or reject the move
      accept = random().metropolis(rosen_f/rosen_r);
      if (accept) {
         // If accepted, keep current positions.
         // Increment counter for accepted moves.
         incrementNAccept();
      } else {
         // If rejected, restore original atom positions
         atomId = beginId;
         for (i = 0; i < nRegrow_; ++i) {
            atom0Ptr = &(molPtr->atom(atomId));
            atom0Ptr->position() = oldPos_[i];
            #ifndef SIMP_NOPAIR
            system().pairPotential().updateAtomCell(*atom0Ptr);
            #endif
            atomId += sign;
         }
      }

      return accept;
   }
}
