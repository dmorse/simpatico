#ifndef LINK_LENGTH_DIST_CPP
#define LINK_LENGTH_DIST_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "LinkLengthDist.h"
#include <mcMd/simulation/Simulation.h>
#include <mcMd/links/LinkMaster.h>
#include <mcMd/species/Species.h>
#include <util/boundary/Boundary.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <mcMd/util/FileMaster.h>

#include <util/global.h>

namespace McMd
{

   using namespace Util;

   /// Constructor.
   LinkLengthDist::LinkLengthDist(System& system) 
    : SystemDiagnostic<System>(system)
   {  setClassName("LinkLengthDist"); }

   /// Read parameters from file, and allocate data array.
   void LinkLengthDist::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);
      read<int>(in, "speciesId", speciesId_);
      if (speciesId_ < 0) {
         UTIL_THROW("Negative speciesId");
      }
      readParamComposite(in, accumulator_);
   }

   /*
   * Clear accumulator.
   */
   void LinkLengthDist::setup() 
   {  accumulator_.clear(); }
 
   /// Add particle pairs to RDF histogram.
   void LinkLengthDist::sample(long iStep) 
   {
      if (isAtInterval(iStep))  {
	 Link     link;
	 int      nLinks0, idLink;	 
         double   lsq, l;
        
         // Go through all the links.
	 nLinks0 = system().linkMaster().nLink();
         for (idLink=0; idLink < nLinks0; idLink++) {         
            link = system().linkMaster().link(idLink);
	    lsq = system().boundary().distanceSq(link.atom0().position(), link.atom1().position());
            l = sqrt(lsq);
            accumulator_.sample(l);         
	 }      
         
      }
   }  


   /// Output results to file after simulation is completed.
   void LinkLengthDist::output() 
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
