/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "InterIntraLink.h"
#include <mcMd/simulation/Simulation.h>
#include <mcMd/links/LinkMaster.h>        
#include <util/misc/FileMaster.h>  
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>


namespace McMd
{

   using namespace Util;

   /// Constructor.
   InterIntraLink::InterIntraLink(System& system) 
    : SystemAnalyzer<System>(system)
   {  setClassName("InterIntraLink"); }

   /// Read parameters from file, and allocate data array.
   void InterIntraLink::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);
      //read<int>(in,"nSamplePerBlock", nSamplePerBlock_);

      nSamplePerBlock_ = 0;
      accumulatorInter_.setNSamplePerBlock(nSamplePerBlock_);
      accumulatorIntra_.setNSamplePerBlock(nSamplePerBlock_);

      // If nSamplePerBlock != 0, open an output file for block averages.
      if (nSamplePerBlock_ != 0) {
         fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      }
      
   }

   /*
   * Clear accumulator.
   */
   void InterIntraLink::setup() 
   {  accumulatorInter_.clear(); 
      accumulatorIntra_.clear();
   }
 
   /// Evaluate inter and intramolecule number of links.
   void InterIntraLink::sample(long iStep) 
   { 
      int       interLink, intraLink;
      int       nLink, iLink, iMol0, iMol1;
      Link*     linkPtr;
      Atom      *atom0Ptr, *atom1Ptr;
      Molecule  *mol0Ptr, *mol1Ptr;
     
      if (isAtInterval(iStep))  {
	
	interLink = 0;
        intraLink = 0;
	
        nLink = system().linkMaster().nLink();
        for (iLink = 0; iLink < nLink ; ++iLink) { 	
	  linkPtr = &(system().linkMaster().link(iLink)); 
	  atom0Ptr = &(linkPtr->atom0());
	  atom1Ptr = &(linkPtr->atom1());
	  mol0Ptr = &atom0Ptr->molecule();
	  mol1Ptr = &atom1Ptr->molecule();
	  iMol0 = system().moleculeId(*mol0Ptr);
	  iMol1 = system().moleculeId(*mol1Ptr);
	  if (iMol0==iMol1) {
	    intraLink+=1;	    
	  }
	  else {
	    interLink+=1;
	  }
        }         
        accumulatorInter_.sample(interLink, outputFile_);
	accumulatorIntra_.sample(intraLink, outputFile_);
      } 
   }

   /*
   * Output results to file after simulation is completed.
   */
   void InterIntraLink::output() 
   { 
      // If outputFile_ was used to write block averages, close it.
      if (nSamplePerBlock_ != 0) {
         outputFile_.close();
      }

      // Write parameters to file
      fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
      ParamComposite::writeParam(outputFile_); 
      outputFile_.close();

      // Write average to file
      fileMaster().openOutputFile(outputFileName(".ave"), outputFile_);
      outputFile_ << "Average number of intermolecular links: " << std::endl;
      accumulatorInter_.output(outputFile_);
      outputFile_ << std::endl;
      outputFile_ << "Average number of intramolecular links: " << std::endl;      
      accumulatorIntra_.output(outputFile_);
      outputFile_.close();
   }

}
