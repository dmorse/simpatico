/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "HybridNphMdMove.h"
#include <mcMd/mcSimulation/McSystem.h>
#include <mcMd/mdIntegrators/MdIntegrator.h>
#include <mcMd/mdSimulation/MdSystem.h>
#include <mcMd/simulation/Simulation.h>
#include <util/boundary/OrthorhombicBoundary.h>
#include <util/ensembles/BoundaryEnsemble.h>
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
   HybridNphMdMove::HybridNphMdMove(McSystem& system) :
      SystemMove(system),
      mdSystemPtr_(0),
      nStep_(0),
      nphIntegratorPtr_(0),
      barostatMass_(0.0),
      mode_()
   {
      setClassName("HybridNphMdMove");
      mdSystemPtr_ = new MdSystem(system);
      oldPositions_.allocate(simulation().atomCapacity());
   }

   /*
   * Destructor.
   */
   HybridNphMdMove::~HybridNphMdMove()
   {
      if (mdSystemPtr_) {
         delete mdSystemPtr_;
      }
   }

   /*
   * Read parameter maxDisp
   */
   void HybridNphMdMove::readParameters(std::istream& in)
   {
      readProbability(in);
      read<int>(in, "nStep", nStep_);
      readParamComposite(in, *mdSystemPtr_);
      
      nphIntegratorPtr_ = dynamic_cast<NphIntegrator*>(&mdSystemPtr_->mdIntegrator());
      
      barostatMass_ = nphIntegratorPtr_->barostatMass();
      mode_ = nphIntegratorPtr_->mode();
   }

   /*
   * Load internal state from an archive.
   */
   void HybridNphMdMove::loadParameters(Serializable::IArchive &ar)
   {
      McMove::loadParameters(ar);
      loadParameter<int>(ar, "nStep", nStep_);
      loadParamComposite(ar, *mdSystemPtr_);
      nphIntegratorPtr_ = dynamic_cast<NphIntegrator*>(&mdSystemPtr_->mdIntegrator());
      
      barostatMass_ = nphIntegratorPtr_->barostatMass();
      mode_ = nphIntegratorPtr_->mode();
   }

   /*
   * Save internal state to an archive.
   */
   void HybridNphMdMove::save(Serializable::OArchive &ar)
   {
      McMove::save(ar);
      ar << nStep_;
      mdSystemPtr_->saveParameters(ar);
   }


   /*
   * Generate, attempt and accept or reject a Hybrid MD/MC move.
   */
   bool HybridNphMdMove::move()
   {
      System::MoleculeIterator molIter;
      Molecule::AtomIterator   atomIter;
      double oldEnergy, newEnergy;
      int    iSpec;
      int    nSpec = simulation().nSpecies();

      bool   accept;

      if (nphIntegratorPtr_ == NULL) {
         UTIL_THROW("null integrator pointer");
      }
      // Increment counter for attempted moves
      incrementNAttempt();
      
      // Store old boundary lengths.
      Vector oldLengths = system().boundary().lengths();

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
      nphIntegratorPtr_->setup();
      
      // generate integrator variables from a Gaussian distribution
      Random& random = simulation().random();
      
      double temp = system().energyEnsemble().temperature();
       
      double volume = system().boundary().volume();
      
      if (mode_ == Cubic) {
         // one degree of freedom
	 // barostat_energy = 1/2 (1/W) eta_x^2
         double sigma = sqrt(temp/barostatMass_);
         nphIntegratorPtr_->setEta(0, sigma*random.gaussian());
      } else if (mode_ == Tetragonal) {
         // two degrees of freedom
         // barostat_energy = 1/2 (1/W) eta_x^2 + 1/2 (1/(2W)) eta_y^2
         double sigma1 = sqrt(temp/barostatMass_);
         nphIntegratorPtr_->setEta(0, sigma1*random.gaussian());
         double sigma2 = sqrt(temp/barostatMass_/2.0);
         nphIntegratorPtr_->setEta(1, sigma2*random.gaussian());
      } else if (mode_ == Orthorhombic) { 
         // three degrees of freedom 
         // barostat_energy = 1/2 (1/W) (eta_x^2 + eta_y^2 + eta_z^2)
         double sigma = sqrt(temp/barostatMass_);
         nphIntegratorPtr_->setEta(0, sigma*random.gaussian());
         nphIntegratorPtr_->setEta(1, sigma*random.gaussian());
         nphIntegratorPtr_->setEta(2, sigma*random.gaussian());
      }

      // Store old energy
      oldEnergy  = mdSystemPtr_->potentialEnergy();
      oldEnergy += mdSystemPtr_->kineticEnergy();
      oldEnergy += system().boundaryEnsemble().pressure()*volume;
      oldEnergy += nphIntegratorPtr_->barostatEnergy();

      // Run a short MD simulation
      for (int iStep = 0; iStep < nStep_; ++iStep) {
         nphIntegratorPtr_->step();
      }
      
      volume = system().boundary().volume();

      // Calculate new energy
      newEnergy  = mdSystemPtr_->potentialEnergy();
      newEnergy += mdSystemPtr_->kineticEnergy();
      newEnergy += system().boundaryEnsemble().pressure()*volume;
      newEnergy += nphIntegratorPtr_->barostatEnergy();

      // Decide whether to accept or reject
      accept = random.metropolis( boltzmann(newEnergy-oldEnergy) );

      // Accept move
      if (accept) {
         
         #ifndef SIMP_NOPAIR
         // Rebuild the McSystem cellList using the new positions.
         system().pairPotential().buildCellList();
         #endif

         // Increment counter for the number of accepted moves.
         incrementNAccept();

      } else {
         
         // Restore old boundary lengths
         system().boundary().setOrthorhombic(oldLengths);

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
