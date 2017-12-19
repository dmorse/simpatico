#ifndef SS_CHAIN_DIST_CPP
#define SS_CHAIN_DIST_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SSChainDist.h"
#include <mcMd/simulation/Simulation.h>
#include <mcMd/links/LinkMaster.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <simp/species/Species.h>
#include <simp/boundary/Boundary.h>
#include <util/misc/FileMaster.h>

#include <util/global.h>

namespace McMd
{

   using namespace Util;
   using namespace Simp;

   /// Constructor.
   SSChainDist::SSChainDist(System& system) 
    : SystemAnalyzer<System>(system)
   {  setClassName("SSChainDist"); }

   /// Read parameters from file.
   /// Suggested values for the histogram:
   /// min   -1.0
   /// max    chainlength                        
   /// nBin   chainlength+1
   
   void SSChainDist::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);
      readParamComposite(in, accumulator_);
   }

   
   /// Clear accumulator.   
   void SSChainDist::setup() 
   {  accumulator_.clear(); }
 
   /// Add atoms attached to links to SSChainDist histogram.
   void SSChainDist::sample(long iStep) 
   {
      if (isAtInterval(iStep))  {
	 Link     link;
	 int      nLinks0, idLink;
	 int      iAtom0, iAtom1;
        
         // Go through all the links.
	 nLinks0 = system().linkMaster().nLink();
         for (idLink=0; idLink < nLinks0; idLink++) {         
            link = system().linkMaster().link(idLink);
	    iAtom0 = link.atom0().indexInMolecule();
	    iAtom1 = link.atom1().indexInMolecule();
            accumulator_.sample((double)iAtom0); 
	    accumulator_.sample((double)iAtom1);
	 }      
         
      }
   }  


   /// Output results to file after simulation is completed.
   void SSChainDist::output() 
   {  

      // Echo parameters to a log file
      fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
      writeParam(outputFile_); 
      outputFile_.close();

      // Output statistical analysis to separate data file
      fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      accumulator_.output(outputFile_); 
      outputFile_.close();

   }

}
#endif
