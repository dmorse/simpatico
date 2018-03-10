#ifndef GC_SLIPLINKMOVE_CPP
#define GC_SLIPLINKMOVE_CPP

/*
* MolMcD - Monte Carlo and Molecular Dynamics Simulator for Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "GcSliplinkMove.h"
#include <util/misc/FileMaster.h>
#include <mcMd/simulation/Simulation.h>
#include <mcMd/links/LinkMaster.h>
#include <mcMd/mcSimulation/McSystem.h>
#include <mcMd/potentials/bond/BondPotential.h>
#include <mcMd/potentials/pair/McPairPotential.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <simp/boundary/Boundary.h>
#include <util/space/Vector.h>
#include <util/global.h>

#include <cstdlib>

#define MCMD_LINK_CREATE
#define MCMD_LINK_SHUFFLE

namespace McMd
{

   using namespace Util;
   using namespace Simp;

   // Constructor
   GcSliplinkMove::GcSliplinkMove(McSystem& system) :
      SystemMove(system),
      cutoff_(0.0),
      cutoffSq_(0.0),
      fCreate_(0.0),
      fNotCreate_(1.0),
      mu_(0.0),
      nTrial_(0),
      speciesId_(0)
   {  setClassName("GcSliplinkMove"); }

   void GcSliplinkMove::readParameters(std::istream& in)
   {
      readProbability(in);
      read<int>(in, "nTrial", nTrial_);
      read<int>(in, "speciesId", speciesId_);
      read<double>(in, "fCreate", fCreate_);
      read<double>(in, "cutoff", cutoff_);
      read<double>(in, "mu", mu_);

      cutoffSq_   = cutoff_*cutoff_;
      fNotCreate_ = 1.0 - fCreate_;
   }

   /// Create or destroy the slip-springs.
   bool GcSliplinkMove::move()
   {
      System::MoleculeIterator molIter;
      double    rsq, oldEnergy, newEnergy;
      double    rnd, norm, energy, sum, pci, pdi, ratio;
      double    cdf[CellList::MaxNeighbor];
      Molecule *mol0Ptr, *mol1Ptr;
      Link     *linkPtr;
      Atom     *atom0Ptr, *atom1Ptr;
      int       i, j, iAtom0, iAtom1, linkId, nLink, nNeighbor, n0, endId;
      int       idneighbors[CellList::MaxNeighbor];
      bool      allowed;

      
      // Loop over attempted Markov moves
      for (i=0; i < nTrial_; ++i) {

         nLink = system().linkMaster().nLink();

         incrementNAttempt();
         if (random().uniform(0.0, 1.0) > fCreate_) { // Shuffle/destroy

            if (nLink > 0) {
   
               // Choose a link at random
               linkId = random().uniformInt(0, nLink);
               linkPtr = &(system().linkMaster().link(linkId));
   
               // Choose an end to be moved: atom0.
               if (random().uniform(0.0, 1.0) < 0.5) {
                  endId = 0;
                  atom0Ptr = &(linkPtr->atom0());
                  atom1Ptr = &(linkPtr->atom1());
               } else {
                  endId = 1;
                  atom0Ptr = &(linkPtr->atom1());
                  atom1Ptr = &(linkPtr->atom0());
               }
   
               // Identify molecules and indices within molecules
               mol0Ptr = &atom0Ptr->molecule();
               iAtom0 = atom0Ptr->indexInMolecule();
               mol1Ptr = &atom1Ptr->molecule();
               iAtom1 = atom1Ptr->indexInMolecule();
   
               rsq = boundary().distanceSq(atom0Ptr->position(),
                                           atom1Ptr->position());
   
               // If atom0 is not a chain end, attempt shuffle
               if (iAtom0 != 0 && iAtom0 != mol0Ptr->nAtom() - 1) {
 
                  #ifdef MCMD_LINK_SHUFFLE 
                  // Calculate old link energy 
                  oldEnergy = system().linkPotential()
                                      .energy(rsq, linkPtr->typeId());
   
                  // Shuffle atom0 index
                  if (random().uniform(0.0, 1.0) > 0.5) {
                     ++iAtom0;
                  } else {
                     --iAtom0;
                  }
                  atom0Ptr = &mol0Ptr->atom(iAtom0);
   
                  // Check if proposed shuffle is allowed
                  allowed = true;
                  if (mol0Ptr == mol1Ptr) {
                     if (std::abs(iAtom0 - iAtom1) < 2) {
                        allowed = false;
                     }
                  }
   
                  if (allowed) {
   
                     // Calculate energy of new configuration
                     rsq = boundary().distanceSq(atom0Ptr->position(),
                                                 atom1Ptr->position());
                     newEnergy = system().linkPotential()
                                         .energy(rsq, linkPtr->typeId());
      
                     // Accept or reject shuffle
                     if (random().metropolis(boltzmann(newEnergy - oldEnergy))) {
                        system().linkMaster().reSetAtom(*linkPtr, *atom0Ptr, endId);
                        //system().linkMaster().reSetAtoms(*linkPtr,
                        //                                 *atom0Ptr, *atom1Ptr);
                        incrementNAccept();	
                     }
   
                  }
                  #endif
   
               } else { // If  atom0 is a chain end
   
                  // Attempt shuffle away from end with 50% probability
                  if (random().uniform(0.0, 1.0) > 0.5) {
   
                     #ifdef MCMD_LINK_SHUFFLE 
                     // Calculate energy of old configuration.
                     oldEnergy = system().linkPotential()
                                         .energy(rsq, linkPtr->typeId());	
   	
                     // Move the atom away from the end
                     if (iAtom0 == 0) {
                        ++iAtom0;
                     } else {
                        assert(iAtom0 == mol0Ptr->nAtom()-1);
                        --iAtom0;
                     }
                     atom0Ptr = &mol0Ptr->atom(iAtom0);	
   
                     // Check if proposed shuffle is allowed
                     allowed = true;
                     if (mol0Ptr == mol1Ptr) {
                        if (std::abs(iAtom0 - iAtom1) < 2) {
                           allowed = false;
                        }
                     }
   
                     if (allowed) {
   
                        // Calculate new link energy 
                        rsq = boundary().distanceSq(atom0Ptr->position(),
                                                    atom1Ptr->position());
                        newEnergy = system().linkPotential()
                                            .energy(rsq, linkPtr->typeId());
      
                        // Accept or reject shuffle
                        if (random().metropolis(boltzmann(newEnergy-oldEnergy))) {
                           system().linkMaster().reSetAtom(*linkPtr, *atom0Ptr, endId);
                           //system().linkMaster()
                           //        .reSetAtoms(*linkPtr, *atom0Ptr, *atom1Ptr);
                           incrementNAccept();
                        }

                     }
                     #endif
    
                  } else { // Attempt destroy with 50% probability
  
                     #ifdef MCMD_LINK_CREATE
                     // If link length is within create/destroy cutoff 
                     if (rsq <= cutoffSq_) {
   
                        // Get array of neighbors
                        system().pairPotential().cellList()
                                .getNeighbors(atom0Ptr->position(), neighbors_);
                        nNeighbor = neighbors_.size();
      
                        // Loop over neighboring atoms
                        n0 = 0;
                        sum = 0;
                        for (j = 0; j < nNeighbor; ++j) {
                           atom1Ptr = neighbors_[j];
      
                           // If atoms are not the same
                           if (atom0Ptr != atom1Ptr) {

                              // If not a masked pair
                              if (!atom0Ptr->mask().isMasked(*atom1Ptr)) {
      
                                 rsq = boundary().distanceSq(atom0Ptr->position(),
                                                              atom1Ptr->position());

                                 // If r < cutoff, compute Boltzmann weight,
                                 // and add to cumulative distribution cdf
                                 if (rsq <= cutoffSq_) {
                                    energy = system().linkPotential()
                                                     .energy(rsq, linkPtr->typeId());
                                    sum = sum + boltzmann(energy);
                                    cdf[n0] = sum;
                                    idneighbors[n0] = j;
                                    n0++;
                                 }

                              }
                           }
                        }
      
                        // Accept or reject destruction
                        pci = 2.0 * system().nMolecule(speciesId_)
                                 * boltzmann(-mu_) * cdf[n0-1] / fCreate_;
                        pdi = 4.0 * nLink / fNotCreate_;
                        ratio = pdi / pci;
                        if (random().metropolis(ratio)) {
                           system().linkMaster().removeLink(linkId);
                           incrementNAccept();
                        }
      
                     } // end if dRsq < cutoffSq
                     #endif

                  } // end attempt destroy
   
               } // end if atom0 at end

            } // end if nLink > 0

         } else { // Attempt slip-link creation with probability fCreate_

            #ifdef MCMD_LINK_CREATE
            // Choose a molecule at random
            mol0Ptr  = &(system().randomMolecule(speciesId_));
  
            // Choose either molecule end at random
            if (random().uniform(0.0, 1.0) > 0.5) {
               iAtom0 = 0;
            } else {
               iAtom0 = mol0Ptr->nAtom() - 1;
            }
            atom0Ptr = &mol0Ptr->atom(iAtom0);
   
            // Get array of neighbors
            system().pairPotential().cellList()
                    .getNeighbors(atom0Ptr->position(), neighbors_);
            nNeighbor = neighbors_.size();
   
            // Loop over neighboring atoms
            n0 = 0;
            sum = 0;
            for (j = 0; j < nNeighbor; ++j) {
               atom1Ptr = neighbors_[j];
   
               // If atoms are not the same
               if (atom0Ptr != atom1Ptr){

                  // If not a masked pair
                  if (!atom0Ptr->mask().isMasked(*atom1Ptr)) {

                     rsq = boundary().distanceSq(atom0Ptr->position(),
                                                           atom1Ptr->position());

                     // If r < cutoff, compute Boltzmann weight for partner,
                     // and add to cumulative distribution cdf.
                     if (rsq <= cutoffSq_) {
                        energy = system().linkPotential().energy(rsq, 0);
                        sum = sum + boltzmann(energy);
                        cdf[n0] = sum;
                        idneighbors[n0] = j;
                        n0++;
                     }

                  }
               }
            }
   
            // If at least 1 candidate has been found.
            if (n0 > 0) {
   
               // Choose a partner with probability cdf[j]/cdf[n0-1]
               j = 0;
               rnd = random().uniform(0.0, 1.0);
               norm = 1.0/cdf[n0-1];
               while (rnd > cdf[j]*norm ){
                  j = j + 1;
               }
               atom1Ptr = neighbors_[idneighbors[j]];
   
               // Accept or reject creation
               pci = 2.0 * system().nMolecule(speciesId_)
                        * boltzmann(-mu_) * cdf[n0-1] / fCreate_;
               pdi = 4.0 * (nLink + 1.0) / fNotCreate_;
               ratio = pci / pdi;
               if (random().metropolis(ratio)) {
                  if (random().uniform(0.0, 1.0) > 0.5){
                     system().linkMaster().addLink(*atom0Ptr, *atom1Ptr, 0);
                  } else {
                     system().linkMaster().addLink(*atom1Ptr, *atom0Ptr, 0);
                  }
                  incrementNAccept();
               }
   
            } // end if (n0 > 0)
            #endif

         } // end if (shuffle/destroy) else (create)

      } // foreach trial

      return true;
   }

}

#ifdef MCMD_LINK_SHUFFLE
#undef MCMD_LINK_SHUFFLE
#endif

#ifdef MCMD_LINK_CREATE
#undef MCMD_LINK_CREATE
#endif

#endif
