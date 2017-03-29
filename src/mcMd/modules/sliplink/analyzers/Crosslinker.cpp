#ifndef CROSSLINKER_CPP
#define CROSSLINKER_CPP

/*
* MolMcD - Monte Carlo and Molecular Dynamics Simulator for Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Crosslinker.h"
#include <util/format/Dbl.h>
#include <util/format/Int.h>
#include <util/misc/FileMaster.h>
#include <mcMd/simulation/Simulation.h>
#include <mcMd/links/LinkMaster.h>

#include <sstream>

namespace McMd
{
  
   using namespace Util;

   // Constructor
   Crosslinker::Crosslinker(System& system) :
      SystemAnalyzer<System>(system),
      nSample_(0),
      cutoff_(0),
      probability_(0)
   {  setClassName("CrossLinker"); }

   void Crosslinker::readParameters(std::istream& in)
   {
      readInterval(in);
      readOutputFileName(in);
      read<double>(in, "cutoff", cutoff_);
      read<double>(in, "probability", probability_);

   }

   void Crosslinker::setup() 
   {
     // Initialize the private CellList
     cellList_.setAtomCapacity(system().nAtom());
     cellList_.setup(system().boundary(), cutoff_);

   }
 
   // Create links and print.
   void Crosslinker::sample(long iStep) 
   {
      if (isAtInterval(iStep)) {

         // Clear the cellList
         cellList_.clear();
         // Add every atom in this System to the CellList
         System::MoleculeIterator molIter;
         Atom*                     atomPtr;
         for (int iSpec=0; iSpec < system().simulation().nSpecies(); ++iSpec) {
            for (system().begin(iSpec, molIter); molIter.notEnd(); ++molIter) {
               for (int ia=0; ia < molIter->nAtom(); ++ia) {
                  atomPtr = &molIter->atom(ia);
                  system().boundary().shift(atomPtr->position());
                  cellList_.addAtom(*atomPtr);
               }
            }
         }

         //Use the cell list to find neighbours and create links
         CellList::NeighborArray cellNeighbor;
         Vector  iPos, jPos;
         Atom   *iAtomPtr, *jAtomPtr;
         double  dRSq, cutoffSq=cutoff_*cutoff_;
         int     nCellNeighbor, nCellAtom, totCells;
         int     ic, ip, iAtomId, jp, jAtomId;
         // Loop over cells containing primary atom. ic = cell index
         totCells = cellList_.totCells();
         for (ic = 0; ic < totCells; ++ic) {
           // Get Array cellNeighbor of Ids of neighbor atoms for cell ic.
           // Elements 0,..., nCellAtom - 1 contain Ids for atoms in cell ic.
           // Elements nCellAtom,..., nCellNeighbor-1 are from neighboring cells.
           cellList_.getCellNeighbors(ic, cellNeighbor, nCellAtom);
           nCellNeighbor = cellNeighbor.size();
           // Loop over atoms in cell ic
           for (ip = 0; ip < nCellAtom; ++ip) {
              iAtomPtr = cellNeighbor[ip]; 
              iPos     = iAtomPtr->position();
              iAtomId  = iAtomPtr->id();
              // Loop over atoms in all neighboring cells, including cell ic.
              for (jp = 0; jp < nCellNeighbor; ++jp) {
                 jAtomPtr = cellNeighbor[jp]; 
                 jPos     = jAtomPtr->position();
                 jAtomId  = jAtomPtr->id();
                 // Avoid double counting: only count pairs with jAtomId > iAtomId
                 if ( jAtomId > iAtomId ) {
                   // Exclude bonded pairs
                   if (!iAtomPtr->mask().isMasked(*jAtomPtr))  {
                      // Calculate distance between atoms i and j
                      dRSq = system().boundary().distanceSq(iPos, jPos);
                      if (dRSq < cutoffSq) {
                         if(system().simulation().random().uniform(0.0, 1.0) < probability_){ 
                          //create a link between i and j atoms
                           system().linkMaster().addLink(*iAtomPtr, *jAtomPtr, 0);
                         }
                      }
                   }
                 } // end if jAtomId > iAtomId
              } // end for jp (j atom)
           } // end for ip (i atom)
         } // end for ic (i cell)

         // Construct a string representation of integer nSample
         std::stringstream ss;
         std::string       nSampleString;
         ss << nSample_;
         nSampleString = ss.str();

         // Construct new fileName: outputFileName + char(nSample)
         std::string filename;
         filename  = outputFileName();
         filename += nSampleString; 

         // Open output file, write data, and close file
         fileMaster().openOutputFile(filename, outputFile_);
         system().writeConfig(outputFile_);
         outputFile_.close();
         nSample_++;

         // Clear the LinkMaster
         system().linkMaster().clear();

      } // end isAtInterval
   }



}

#endif

