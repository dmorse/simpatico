#ifndef MCMD_MD_MOVE_CPP
#define MCMD_MD_MOVE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "MdMove.h"
#include <mcMd/mcSimulation/McSystem.h>
#include <mcMd/mdSimulation/MdSystem.h>
#include <mcMd/mdIntegrators/MdIntegrator.h>
#ifndef INTER_NOPAIR
#include <mcMd/potentials/pair/MdPairPotential.h>
#include <mcMd/potentials/pair/McPairPotential.h>
#endif

namespace McMd
{

   using namespace Util;

   /*
   * Constructor
   */
   MdMove::MdMove(McSystem& system) :
      SystemMove(system),
      mdSystemPtr_(0),
      nStep_(0)
   {
      mdSystemPtr_ = new MdSystem(system);
   }

   /*
   * Destructor.
   */
   MdMove::~MdMove()
   {
      if (mdSystemPtr_) {
         delete mdSystemPtr_;
      }
   }

   /*
   * Read parameter maxDisp
   */
   void MdMove::readParam(std::istream& in)
   {
      readProbability(in);
      read<int>(in, "nStep", nStep_);
      readParamComposite(in, *mdSystemPtr_);
   }

   /*
   * Generate, attempt and accept or reject a Hybrid MD/MC move.
   */
   bool MdMove::move()
   {
      // Increment counter for attempted moves
      incrementNAttempt();

      // Initialize MdSystem
      #ifndef INTER_NOPAIR
      mdSystemPtr_->pairPotential().buildPairList();
      #endif
      mdSystemPtr_->calculateForces();
      mdSystemPtr_->setBoltzmannVelocities(energyEnsemble().temperature());

      // Run a short MD simulation
      for (int iStep = 0; iStep < nStep_; ++iStep) {
         mdSystemPtr_->mdIntegrator().step();
      }

      #ifndef INTER_NOPAIR
      // Rebuild the McSystem cellList using the new positions.
      system().pairPotential().buildCellList();
      #endif

      // Accept all moves
      incrementNAccept();
      return true;

   }

}
#endif
