#ifndef ATOM_DISPLACE_MOVE_CPP
#define ATOM_DISPLACE_MOVE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "AtomDisplaceMove.h"
#include <mcMd/mcSimulation/McSystem.h>
#ifndef MCMD_NOPAIR
#include <mcMd/potentials/pair/McPairPotential.h>
#endif
#include <mcMd/boundary/Boundary.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <util/space/Vector.h>
#include <util/space/Dimension.h>
#include <util/global.h>

namespace McMd
{

   using namespace Util;

   /* 
   * Constructor
   */
   AtomDisplaceMove::AtomDisplaceMove(McSystem& system) 
    : SystemMove(system),
      delta_(0.0),
      speciesId_(-1)
   {}
   
   /* 
   * Read speciesId and delta.
   */
   void AtomDisplaceMove::readParam(std::istream& in) 
   {
      readProbability(in);
      read<int>(in, "speciesId", speciesId_);
      read<double>(in, "delta", delta_);
   }
   
   /* 
   * Generate, attempt and accept or reject a move.
   */
   bool AtomDisplaceMove::move() 
   { 
      Vector    oldPos;
      double    newEnergy, oldEnergy;
      Molecule* molPtr;
      Atom*     atomPtr;
      int       iAtom;

      incrementNAttempt();

      // Choose a molecule and atom at random
      molPtr  = &(system().randomMolecule(speciesId_));
      iAtom   = random().getInteger(0, molPtr->nAtom());
      atomPtr = &molPtr->atom(iAtom);

      // Calculate current pair energy and store old position
      oldPos    = atomPtr->position();
      oldEnergy = system().atomPotentialEnergy(*atomPtr);

      for (int j = 0; j < Dimension; ++j) {
         atomPtr->position()[j] += random().getFloat(-delta_, delta_);
      }
      boundary().shift(atomPtr->position());
      newEnergy = system().atomPotentialEnergy(*atomPtr);

      // Decide whether to accept forward move
      bool accept = random().metropolis(boltzmann(newEnergy - oldEnergy));

      if (accept) {
  
         #ifndef MCMD_NOPAIR
         system().pairPotential().updateAtomCell(*atomPtr);
         #endif
         incrementNAccept();

      } else {
   
         atomPtr->position() = oldPos;
   
      }
      return accept;
   }

}      
#endif
