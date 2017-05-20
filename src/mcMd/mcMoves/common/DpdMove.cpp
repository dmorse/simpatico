/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "DpdMove.h"
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
   DpdMove::DpdMove(McSystem& system) :
      SystemMove(system),
      mdSystemPtr_(0),
      nStep_(0)
   {
      setClassName("DpdMove");
      mdSystemPtr_ = new MdSystem(system);
   }

   /*
   * Destructor.
   */
   DpdMove::~DpdMove()
   {
      if (mdSystemPtr_) {
         delete mdSystemPtr_;
      }
   }

   /*
   * Read probability, nStep, and MdSystem parameters.
   */
   void DpdMove::readParameters(std::istream& in)
   {
      readProbability(in);
      read<int>(in, "nStep", nStep_);
      readParamComposite(in, *mdSystemPtr_);
   }

   /*
   * Load internal state from an archive.
   */
   void DpdMove::loadParameters(Serializable::IArchive &ar)
   {
      McMove::loadParameters(ar);
      loadParameter<int>(ar, "nStep", nStep_);
      loadParamComposite(ar, *mdSystemPtr_);
   }

   /*
   * Save internal state to an archive.
   */
   void DpdMove::save(Serializable::OArchive &ar)
   {
      McMove::save(ar);
      ar << nStep_;
      mdSystemPtr_->saveParameters(ar);
   }

   /*
   * Setup before simulation run.
   */
   void DpdMove::setup()
   {  mdSystemPtr_->mdIntegrator().setup(); }

   /*
   * Generate, attempt and accept or reject a Hybrid MD/MC move.
   */
   bool DpdMove::move()
   {
      // Increment counter for attempted moves
      incrementNAttempt();
      
      // Calculate conservative forces
      #ifndef SIMP_NOPAIR
      mdSystemPtr_->pairPotential().buildPairList();
      #endif
      mdSystemPtr_->calculateForces();

      // Run a short DPD MD simulation
      for (int iStep = 0; iStep < nStep_; ++iStep) {
         mdSystemPtr_->mdIntegrator().step();
      }

      #ifndef SIMP_NOPAIR
      // Rebuild the McSystem cellList using the new positions.
      system().pairPotential().buildCellList();
      #endif

      // Increment counter for the number of accepted moves.
      incrementNAccept();
      return true;

   }

}
