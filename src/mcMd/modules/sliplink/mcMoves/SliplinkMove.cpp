#ifndef SLIPLINKMOVE_CPP
#define SLIPLINKMOVE_CPP

/*
* MolMcD - Monte Carlo and Molecular Dynamics Simulator for Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SliplinkMove.h"
#include <util/misc/FileMaster.h>
#include <mcMd/simulation/Simulation.h>
#include <mcMd/links/LinkMaster.h>
#include <mcMd/mcSimulation/McSystem.h>
#include <mcMd/potentials/pair/McPairPotential.h>
#include <mcMd/potentials/bond/BondPotential.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <simp/boundary/Boundary.h>
#include <util/space/Vector.h>
#include <util/global.h>
#include <cmath>

#define MCMD_LINK_TRANSFER

namespace McMd
{

   using namespace Util;
   using namespace Simp;

   // Constructor
   SliplinkMove::SliplinkMove(McSystem& system) :
      SystemMove(system),
      cutoff_(0),
      speciesId_(0)
   { setClassName("SliplinkMove"); }

   void SliplinkMove::readParameters(std::istream& in)
   {
      readProbability(in);
      read<double>(in, "cutoff", cutoff_);
      read<int>(in, "speciesId", speciesId_);
   }


   // Create or destroy the slip-springs.
   bool SliplinkMove::move()
   {
      System::MoleculeIterator molIter;
      Atom                    *atom0Ptr, *atom1Ptr;
      Molecule                *mol0Ptr, *mol1Ptr;
      int                      j, iMolecule0, iMolecule1, n0, nLinks0, endId;
      int                      idLink, iAtom0,iAtom1, iAtom, nNeighbor,id0, id1;
      double                   rsq, oldEnergy, newEnergy;
      double                   dRSq, mindRSq=cutoff_*cutoff_, rnd, norm;
      Link*                    linkPtr;
      static const int         maxNeighbor = CellList::MaxNeighbor;
      double                   cdf[maxNeighbor], energy, sum;
      int                      idneighbors[maxNeighbor];


      // Go through all links.
      nLinks0 = system().linkMaster().nLink();
      for (idLink = nLinks0-1; idLink >= 0 ; idLink--) {
         incrementNAttempt();
         linkPtr = &(system().linkMaster().link(idLink));

         // Choose an end to be moved: atom0.
         endId = 0;	 
         atom0Ptr = &(linkPtr->atom0());
         atom1Ptr = &(linkPtr->atom1());
         if (random().uniform(0.0, 1.0) > 0.5) {
	    endId=1;	   
            atom0Ptr = &(linkPtr->atom1());
            atom1Ptr = &(linkPtr->atom0());
         }
         mol0Ptr = &atom0Ptr->molecule();
         iAtom0 = atom0Ptr->indexInMolecule();
         iMolecule0 = system().moleculeId(*mol0Ptr);
         mol1Ptr = &atom1Ptr->molecule();
         iMolecule1 = system().moleculeId(*mol1Ptr);
         iAtom1 = atom1Ptr->indexInMolecule();

         // If the atom is not a chain end
         if (iAtom0 != 0 && iAtom0 != mol0Ptr->nAtom() - 1) {

            // Calculate energy of old configuration.
            rsq = boundary().distanceSq(atom0Ptr->position(), atom1Ptr->position());
            oldEnergy = system().linkPotential().energy(rsq, linkPtr->typeId());

            // Move the atom
            iAtom = iAtom0;
            iAtom0 = iAtom - 1;
            if(random().uniform(0.0, 1.0) > 0.5) iAtom0 = iAtom + 1;
            // If the atoms belong to the same chain
            if(iMolecule0==iMolecule1){
                if(std::fabs(iAtom0-iAtom1)<2) iAtom0 = iAtom;
            }
            atom0Ptr = &mol0Ptr->atom(iAtom0);

            // Calculate energy of new configuration
            rsq = boundary().distanceSq(atom0Ptr->position(), atom1Ptr->position());
            newEnergy = system().linkPotential().energy(rsq, linkPtr->typeId());
            // Accept new configuration with Metropolis
            if (random().metropolis(boltzmann(newEnergy - oldEnergy))) {
	       system().linkMaster().reSetAtom(*linkPtr, *atom0Ptr, endId);
	       //system().linkMaster().reSetAtoms(*linkPtr, *atom0Ptr, *atom1Ptr);
               //system().linkMaster().addLink(*atom0Ptr, *atom1Ptr, 0);
               //system().linkMaster().removeLink(idLink);
               incrementNAccept();
            }

         } else { // If the atom is a chain end

            rsq = boundary().distanceSq(atom0Ptr->position(), atom1Ptr->position());
 
            // try to move in the chain
             if (random().uniform(0.0, 1.0) > 0.5) {
 
               // Calculate energy of old configuration.
               oldEnergy = system().linkPotential().energy(rsq, linkPtr->typeId());

               // Move the atom
               iAtom = iAtom0;
               iAtom0 = iAtom - 1;
               if (iAtom == 0) iAtom0 = 1;
               // If the atoms belong to the same chain
               if(iMolecule0==iMolecule1){
                 if(std::fabs(iAtom0-iAtom1)<2) iAtom0 = iAtom;
               }
               atom0Ptr = &mol0Ptr->atom(iAtom0);

               // Calculate energy of new configuration
               rsq = boundary().distanceSq(atom0Ptr->position(),
                                           atom1Ptr->position());
               newEnergy = system().linkPotential().energy(rsq, linkPtr->typeId());

               // Accept or reject new configuration.
               if (random().metropolis(boltzmann(newEnergy - oldEnergy))) {
   	          system().linkMaster().reSetAtom(*linkPtr, *atom0Ptr, endId);
		  //system().linkMaster().reSetAtoms(*linkPtr, *atom0Ptr, *atom1Ptr);
                  incrementNAccept();
               }

            } else { // try to renew the slip-link

               #ifdef MCMD_LINK_TRANSFER
               // Try to renew the link if the bond is smaller than the cutoff
               if (rsq <= mindRSq) {

                  // Get array of neighbors
                  system().pairPotential().cellList().getNeighbors(atom0Ptr->position(), neighbors_);
                  nNeighbor = neighbors_.size();

                  iMolecule0 = system().moleculeId(*mol0Ptr);
                  id0 = atom0Ptr->id();
                  n0 = 0;
                  sum = 0;

                  // Loop over neighboring atoms
                  for (j = 0; j < nNeighbor; ++j) {
                     atom1Ptr = neighbors_[j];
                     mol1Ptr = &atom1Ptr->molecule();
                     iMolecule1 = system().moleculeId(*mol1Ptr);
                     id1 = atom1Ptr->id();

                     // Check if atoms are the same
                     if (id0 != id1){

                        // Exclude masked pairs
                        if (!atom0Ptr->mask().isMasked(*atom1Ptr)) {

                           // Identify possible partners and compute weights.
                           dRSq = system().boundary().distanceSq(atom0Ptr->position(),
                                                                 atom1Ptr->position());
                           if (dRSq <= mindRSq) {
                              energy = system().linkPotential().energy(dRSq, 0);
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
                     // renew the slip-link with the selected partner.
                     system().linkMaster().addLink(*atom0Ptr, *atom1Ptr, 0);
                     system().linkMaster().removeLink(idLink);
                     incrementNAccept();
                  }

               } // else if dRsq < mindRSq
               #endif

            } // if (shuffle) else (transfer)

         } // if (not at end) else (at end)

      } // foreach link
      return true;
   }

}
#endif
