#ifndef MCMD_MC_WRITE_RESTART_CPP
#define MCMD_MC_WRITE_RESTART_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "McWriteRestart.h"
#include <mcMd/mcSimulation/McSimulation.h>
//#include <mcMd/util/FileMaster.h>
#include <util/util/ioUtil.h>

#include <sstream>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   McWriteRestart::McWriteRestart(McSimulation& simulation) 
    : Diagnostic(),
      filename_(),
      simulationPtr_(&simulation)
   {  setClassName("McWriteRestart"); }

   /*
   * Read interval and outputFileName. 
   */
   void McWriteRestart::readParameters(std::istream& in) 
   {
      readInterval(in);
      read<std::string>(in, "fileName", filename_);
   }

   /*
   * Dump configuration to file
   */
   void McWriteRestart::sample(long iStep) 
   {
      if (isAtInterval(iStep))  {

         // Construct new fileName: outputFileName + toString(nSample)
         //std::string filename;
         //filename  = filename_;
         //filename += toString(nSample_);

         // Open output file, write data, and close file
         simulationPtr_->writeRestart(filename_);

      }
   }
  
}
#endif 
