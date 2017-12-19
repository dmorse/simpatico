#ifndef SLIPLINKER_CPP
#define SLIPLINKER_CPP

/*
* MolMcD - Monte Carlo and Molecular Dynamics Simulator for Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Sliplinker.h"
#include <util/misc/FileMaster.h>
#include <mcMd/simulation/Simulation.h>
#include <mcMd/links/LinkMaster.h>
#include <mcMd/mcSimulation/McSystem.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <simp/boundary/Boundary.h>
#include <util/space/Vector.h>
#include <util/space/Dimension.h>
#include <util/global.h>


namespace McMd
{
  
   using namespace Util;
   using namespace Simp;

   // Constructor
   Sliplinker::Sliplinker(McSystem& system) :
      SystemMove(system),
      cutoff_(0),
      mu_(0),
      speciesId_(0)
   {  setClassName("Sliplinker"); }

   void Sliplinker::readParameters(std::istream& in)
   {
      read<double>(in, "cutoff", cutoff_);
      read<double>(in, "mu", mu_);
      read<int>(in, "speciesId", speciesId_);
   }


   // Create or destroy the slip-springs.
   bool Sliplinker::move() 
   {
     System::MoleculeIterator molIter;
     Molecule::AtomIterator   atomIter;
     Atom                    *atom0Ptr, *atom1Ptr;
     Atom*                    atomPtr;
     Molecule                *mol0Ptr, *mol1Ptr;
     Molecule*                molIPtr;
     double                   prob, dRSq, mindRSq=cutoff_*cutoff_, rnd, norm;
     int                      i, ntrials, j, nNeighbor, idLink, iAtom, id1, id0;
     int                      iAtom0, iAtom1, iMolecule0, iMolecule1, n0;
     Link*                    linkPtr;
     static const int         maxNeighbor = CellList::MaxNeighbor;
     double                   cdf[maxNeighbor], energy, sum;
     int                      idneighbors[maxNeighbor];

     
     ntrials = 2 * system().simulation().atomCapacity();
     for (i=0; i < ntrials; ++i){   
       
       //Choose to create or destroy a link with prob 0.5
       if (random().uniform(0.0, 1.0) > 0.5){         
	  // Try to create a link. 
	  incrementNAttempt();		
	  // Choose a molecule and atom at random
	  molIPtr  = &(system().randomMolecule(speciesId_));
	  iMolecule0 = system().moleculeId(*molIPtr);
	  iAtom   = random().uniformInt(0, molIPtr->nAtom());
	  atomPtr = &molIPtr->atom(iAtom);
	  id0 = atomPtr->id();
	  
	  // Get array of neighbors
	  system().cellList().getNeighbors(atomPtr->position(), neighbors_);
	  nNeighbor = neighbors_.size();
	    
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
	      if (!atomPtr->mask().isMasked(*atom1Ptr)) {
		// Identify possible partners and calculate the cumulative distribution function
		dRSq = system().boundary().distanceSq(atomPtr->position(), atom1Ptr->position());
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
	    prob = 2.0 * (system().linkMaster().nLink() + 1.0);
	    prob = system().simulation().atomCapacity() * boltzmann(-mu_) * cdf[n0-1]/ prob ;
	    //prob = 2.0 * system().nMolecule(speciesId_) * boltzmann(-mu_) * cdf[n0-1]/ prob ;	      
	    if (system().simulation().random().uniform(0.0, 1.0) < prob) {        
	      system().linkMaster().addLink(*atomPtr, *atom1Ptr, 0); 
	      incrementNAccept();	    
	    }      
	  }    
       }
       else {
	  // Try to destroy a link
	  incrementNAttempt();
	  // Choose a link at random
	  if (system().linkMaster().nLink() > 0){
	    idLink = random().uniformInt(0, system().linkMaster().nLink());
	    // Indentify the atoms for this link.
	    linkPtr = &(system().linkMaster().link(idLink));
	    atom0Ptr = &(linkPtr->atom0());
	    mol0Ptr = &atom0Ptr->molecule();
	    iMolecule0 = system().moleculeId(*mol0Ptr);
	    iAtom0 = atom0Ptr->indexInMolecule();	    
	    atom1Ptr = &(linkPtr->atom1());
	    mol1Ptr = &atom1Ptr->molecule();
	    iMolecule1 = system().moleculeId(*mol1Ptr);
	    iAtom1 = atom1Ptr->indexInMolecule();
	    // try to delete the slip-spring	
	    dRSq = system().boundary().distanceSq(atom0Ptr->position(), atom1Ptr->position());
	    // try to delete the link if the bond is smaller than the cutoff
	    if (dRSq <= mindRSq) {    
	      // Get array of neighbors
	      system().cellList().getNeighbors(atom0Ptr->position(), neighbors_);
	      nNeighbor = neighbors_.size();
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
	      //prob = 2.0 * system().nMolecule(speciesId_) * boltzmann(-mu_) * cdf[n0-1];	
	      prob = system().simulation().atomCapacity() * boltzmann(-mu_) * cdf[n0-1];	    
	      prob = 2.0 * system().linkMaster().nLink() / prob; 	            
	      if (system().simulation().random().uniform(0.0, 1.0) < prob) {        
		system().linkMaster().removeLink(idLink); 
		incrementNAccept();  
	      }  
	    }
          }
       }
     }
     return true;
   }

}
#endif
