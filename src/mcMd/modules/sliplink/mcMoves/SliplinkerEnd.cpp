#ifndef SLIPLINKER_END_CPP
#define SLIPLINKER_END_CPP

/*
* MolMcD - Monte Carlo and Molecular Dynamics Simulator for Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SliplinkerEnd.h"
#include <mcMd/simulation/Simulation.h>
#include <mcMd/mcSimulation/McSystem.h>
#include <mcMd/links/LinkMaster.h>
#include <mcMd/potentials/pair/McPairPotential.h>
#include <mcMd/potentials/bond/BondPotential.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <simp/boundary/Boundary.h>
#include <util/misc/FileMaster.h>
#include <util/space/Vector.h>
#include <util/space/Dimension.h>
#include <util/global.h>


namespace McMd
{

   using namespace Util;
   using namespace Simp;

   // Constructor
   SliplinkerEnd::SliplinkerEnd(McSystem& system) :
      SystemMove(system),
      cutoff_(0),
      mu_(0),
      speciesId_(0)
   { setClassName("SliplinkerEnd"); }

   void SliplinkerEnd::readParameters(std::istream& in)
   {
      readProbability(in);
      read<double>(in, "cutoff", cutoff_);
      read<double>(in, "mu", mu_);
      read<int>(in, "speciesId", speciesId_);
   }


   // Create or destroy the slip-springs.
   bool SliplinkerEnd::move()
   {
     System::MoleculeIterator molIter;
     Atom                    *atom0Ptr, *atom1Ptr;
     Molecule*                mol0Ptr;
     //Molecule*                mol1Ptr;
     Molecule*                molIPtr;
     double                   prob, dRSq, mindRSq=cutoff_*cutoff_, rnd, norm;
     int                      i, ntrials, j, nNeighbor, idLink, id0, id1;
     int                      iAtom0, n0;
     // int                   iAtom1;
     Link*                    linkPtr;
     static const int         maxNeighbor = CellList::MaxNeighbor;
     double                   cdf[maxNeighbor], energy, sum;
     int                      idneighbors[maxNeighbor];

     ntrials = 4 * system().nMolecule(speciesId_);
     for (i=0; i < ntrials; ++i){

	//Choose to create or destroy a link with prob 0.5
	if (random().uniform(0.0, 1.0) > 0.5){

	    // Try to create a link.
	    incrementNAttempt();

	    // Choose a molecule at random
	    molIPtr  = &(system().randomMolecule(speciesId_));

	    // Choose a molecule end at random
	    iAtom0 = 0;
	    if (random().uniform(0.0, 1.0) > 0.5) iAtom0 = molIPtr->nAtom() - 1;
	    atom0Ptr = &molIPtr->atom(iAtom0);

	    // Get array of neighbors
	    system().pairPotential().cellList().getNeighbors(atom0Ptr->position(), neighbors_);
	    nNeighbor = neighbors_.size();
	    id0 = atom0Ptr->id();

	    n0 = 0;
	    sum = 0;
	    // Loop over neighboring atoms
	    for (j = 0; j < nNeighbor; ++j) {
	      atom1Ptr = neighbors_[j];
	      // mol1Ptr = &atom1Ptr->molecule();
	      id1 = atom1Ptr->id();
	
              // Check if atoms are the same
              if (id0 != id1){	
		// Exclude masked pairs	
		if (!atom0Ptr->mask().isMasked(*atom1Ptr)) {
		  // Identify possible partners and calculate the cumulative distribution function
		  dRSq = system().boundary().distanceSq(atom0Ptr->position(), atom1Ptr->position());
		  if (dRSq <= mindRSq) {
		    energy = system().linkPotential().energy(dRSq, 0);
		    //energy = 0.5*dRSq;
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

	      // Create a slip-link between the selected atoms with probability = prob
	      prob = 2.0*(system().linkMaster().nLink() + 1.0);
	      prob = 2.0 * system().nMolecule(speciesId_) * boltzmann(-mu_) * cdf[n0-1]/ prob ;
	      if (system().simulation().random().uniform(0.0, 1.0) < prob) {
	        if (random().uniform(0.0, 1.0) > 0.5){
		  system().linkMaster().addLink(*atom0Ptr, *atom1Ptr, 0);
                } else {
		  system().linkMaster().addLink(*atom1Ptr, *atom0Ptr, 0);
                }
		incrementNAccept();	
	      }

	    }

	} else {	

	    // Try to destroy a link
	    incrementNAttempt();

	    if (system().linkMaster().nLink() > 0){

	      // Choose a link at random
	      idLink = random().uniformInt(0, system().linkMaster().nLink());
	      linkPtr = &(system().linkMaster().link(idLink));

              // Choose atom0 and atom1 from link ends at random
	      if (random().uniform(0.0, 1.0) > 0.5){
                 atom0Ptr = &(linkPtr->atom0());
                 atom1Ptr = &(linkPtr->atom1());
              } else {
                 atom0Ptr = &(linkPtr->atom1());
                 atom1Ptr = &(linkPtr->atom0());
              }

              mol0Ptr = &atom0Ptr->molecule();
              iAtom0 = atom0Ptr->indexInMolecule();

              // mol1Ptr = &atom1Ptr->molecule();
              // iAtom1 = atom1Ptr->indexInMolecule();

	      // If atom 0 is at the end of a chain, try to delete the slip-spring
	      if (iAtom0 == 0 || iAtom0 == mol0Ptr->nAtom() - 1){

		dRSq = system().boundary().distanceSq(atom0Ptr->position(), atom1Ptr->position());

		// try to delete the link if the bond is smaller than the cutoff
		if (dRSq <= mindRSq) {
		  // Get array of neighbors
		  system().pairPotential().cellList()
                          .getNeighbors(atom0Ptr->position(), neighbors_);
		  nNeighbor = neighbors_.size();
		  id0 = atom0Ptr->id();
		  n0 = 0;
		  sum = 0;
		  // Loop over neighboring atoms
		  for (j = 0; j < nNeighbor; ++j) {
		    atom1Ptr = neighbors_[j];
		    // mol1Ptr = &atom1Ptr->molecule();
	            id1 = atom1Ptr->id();
	
                    // Check if atoms are the same
                    if (id0 != id1){	      		
		      // Exclude masked pairs	
		      if (!atom0Ptr->mask().isMasked(*atom1Ptr)) {
			// Identify possible partners and calculate the cumulative distribution function
			dRSq = system().boundary().distanceSq(atom0Ptr->position(), atom1Ptr->position());
			if (dRSq <= mindRSq) {
			  energy = system().linkPotential().energy(dRSq, 0);
			  //energy = 0.5*dRSq;
			  sum = sum + boltzmann(energy);
			  cdf[n0] = sum;
			  idneighbors[n0] = j;
			  n0++;
			}
		      }
		    }
		  }
		
		  // Destroy the slip-link between the selected atoms with probability = prob	
		  prob = 2.0 * system().nMolecule(speciesId_) * boltzmann(-mu_) * cdf[n0-1];
		  prob = 2.0*system().linkMaster().nLink() / prob; 	
		  if (system().simulation().random().uniform(0.0, 1.0) < prob) {
		    system().linkMaster().removeLink(idLink);
		    incrementNAccept();
		  }
		}
	      }	
	    }
	}
     }
     return true;
   }

}
#endif
