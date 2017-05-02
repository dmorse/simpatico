/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "HybridMdMove.h"
#include <mcMd/mcSimulation/McSystem.h>
#include <mcMd/mdIntegrators/MdIntegrator.h>
#include <mcMd/mdSimulation/MdSystem.h>
#include <mcMd/simulation/Simulation.h>
#include <mcMd/chemistry/Atom.h>
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
   HybridMdMove::HybridMdMove(McSystem& system) :
      SystemMove(system),
      mdSystemPtr_(0),
      nStep_(0)
   {
      setClassName("HybridMdMove");
      mdSystemPtr_ = new MdSystem(system);
      oldPositions_.allocate(simulation().atomCapacity());
   }

   /*
   * Destructor.
   */
   HybridMdMove::~HybridMdMove()
   {
      if (mdSystemPtr_) {
         delete mdSystemPtr_;
      }
   }

   /*
   * Read parameter maxDisp
   */
   void HybridMdMove::readParameters(std::istream& in)
   {
      readProbability(in);
      read<int>(in, "nStep", nStep_);
      readParamComposite(in, *mdSystemPtr_);
   }

   /*
   * Load internal state from an archive.
   */
   void HybridMdMove::loadParameters(Serializable::IArchive &ar)
   {
      McMove::loadParameters(ar);
      loadParameter<int>(ar, "nStep", nStep_);
      loadParamComposite(ar, *mdSystemPtr_);
   }

   /*
   * Save internal state to an archive.
   */
   void HybridMdMove::save(Serializable::OArchive &ar)
   {
      McMove::save(ar);
      ar << nStep_;
      mdSystemPtr_->saveParameters(ar);
   }

   /*
   * Generate, attempt and accept or reject a Hybrid MD/MC move.
   */
   bool HybridMdMove::move()
   {
      System::MoleculeIterator molIter;
      Molecule::AtomIterator   atomIter;
      double oldEnergy, newEnergy;
      int    iSpec;
      int    nSpec = simulation().nSpecies();
      bool   accept;

      incrementNAttempt();

      // Store old atom positions in oldPositions_ array.
      for (iSpec = 0; iSpec < nSpec; ++iSpec) {
         mdSystemPtr_->begin(iSpec, molIter);
         for ( ; molIter.notEnd(); ++molIter) {
            for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
               oldPositions_[atomIter->id()] = atomIter->position();
            }
         }
      }

      // Initialize MdSystem
      #ifndef SIMP_NOPAIR
      mdSystemPtr_->pairPotential().buildPairList();
      #endif
      mdSystemPtr_->calculateForces();
      mdSystemPtr_->setBoltzmannVelocities(energyEnsemble().temperature());
      mdSystemPtr_->mdIntegrator().setup();

      // Store old energy
      oldEnergy  = mdSystemPtr_->potentialEnergy();
      oldEnergy += mdSystemPtr_->kineticEnergy();

      // Run a short MD simulation
      for (int iStep = 0; iStep < nStep_; ++iStep) {
         mdSystemPtr_->mdIntegrator().step();
      }

      // Calculate new energy
      newEnergy  = mdSystemPtr_->potentialEnergy();
      newEnergy += mdSystemPtr_->kineticEnergy();

      // Decide whether to accept or reject
      accept = random().metropolis( boltzmann(newEnergy-oldEnergy) );

      // Accept move
      if (accept) {

         #ifndef SIMP_NOPAIR
         // Rebuild the McSystem cellList using the new positions.
         system().pairPotential().buildCellList();
         #endif

         // Increment counter for the number of accepted moves.
         incrementNAccept();

      } else {

         // Restore old atom positions
         for (iSpec = 0; iSpec < nSpec; ++iSpec) {
            mdSystemPtr_->begin(iSpec, molIter);
            for ( ; molIter.notEnd(); ++molIter) {
               molIter->begin(atomIter);
               for ( ; atomIter.notEnd(); ++atomIter) {
                  atomIter->position() = oldPositions_[atomIter->id()];
               }
            }
         }

      }

      return accept;

   }

}
