/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "MdMove.h"
#include <mcMd/mcSimulation/McSystem.h>
#include <mcMd/mdSimulation/MdSystem.h>
#include <mcMd/mdIntegrators/MdIntegrator.h>
#ifndef SIMP_NOPAIR
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
      setClassName("MdMove");
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
   void MdMove::readParameters(std::istream& in)
   {
      readProbability(in);
      read<int>(in, "nStep", nStep_);
      readParamComposite(in, *mdSystemPtr_);
   }

   /*
   * Load internal state from an archive.
   */
   void MdMove::loadParameters(Serializable::IArchive &ar)
   {
      McMove::loadParameters(ar);
      loadParameter<int>(ar, "nStep", nStep_);
      loadParamComposite(ar, *mdSystemPtr_);
   }

   /*
   * Save internal state to an archive.
   */
   void MdMove::save(Serializable::OArchive &ar)
   {
      McMove::save(ar);
      ar << nStep_;
      mdSystemPtr_->saveParameters(ar);
   }

   /*
   * Generate, attempt and accept or reject a Hybrid MD/MC move.
   */
   bool MdMove::move()
   {
      incrementNAttempt();

      // Initialize MdSystem
      #ifndef SIMP_NOPAIR
      mdSystemPtr_->pairPotential().buildPairList();
      #endif
      mdSystemPtr_->calculateForces();
      mdSystemPtr_->setBoltzmannVelocities(energyEnsemble().temperature());

      // Run a short MD simulation
      for (int iStep = 0; iStep < nStep_; ++iStep) {
         mdSystemPtr_->mdIntegrator().step();
      }

      #ifndef SIMP_NOPAIR
      // Rebuild the McSystem cellList using the new positions.
      system().pairPotential().buildCellList();
      #endif

      // Accept all moves
      incrementNAccept();
      return true;
   }

}
