#ifndef NLINK_AVERAGE_CPP
#define NLINK_AVERAGE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "NLinkAverage.h"
#include <mcMd/simulation/Simulation.h>
#include <mcMd/links/LinkMaster.h>        
#include <mcMd/util/FileMaster.h>        


namespace McMd
{

   using namespace Util;

   /// Constructor.
   NLinkAverage::NLinkAverage(System& system) 
    : SystemDiagnostic<System>(system)
   {}

   /// Read parameters from file, and allocate data array.
   void NLinkAverage::readParam(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);
      read<int>(in,"nSamplePerBlock", nSamplePerBlock_);

      accumulator_.setNSamplePerBlock(nSamplePerBlock_);

      // If nSamplePerBlock != 0, open an output file for block averages.
      if (accumulator_.nSamplePerBlock()) {
         fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      }
   }

   /*
   * Clear accumulator.
   */
   void NLinkAverage::setup() 
   {  accumulator_.clear(); }
 
   /// Evaluate end-to-end vectors of all chains.
   void NLinkAverage::sample(long iStep) 
   { 
      if (isAtInterval(iStep))  {
         double rn = (double) system().linkMaster().nLink();
         accumulator_.sample(rn, outputFile_);
      } 
   }

   /*
   * Output results to file after simulation is completed.
   */
   void NLinkAverage::output() 
   { 
      // If outputFile_ was used to write block averages, close it.
      if (accumulator_.nSamplePerBlock()) {
         outputFile_.close();
      }

      // Write parameters to file
      fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
      ParamComposite::writeParam(outputFile_); 
      outputFile_.close();

      // Write average to file
      fileMaster().openOutputFile(outputFileName(".ave"), outputFile_);
      accumulator_.output(outputFile_); 
      outputFile_.close();
   }

}
#endif 
