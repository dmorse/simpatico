#ifndef MCMD_MD_WRITE_RESTART_CPP
#define MCMD_MD_WRITE_RESTART_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "MdWriteRestart.h"
#include <mcMd/mdSimulation/MdSimulation.h>
//#include <mcMd/misc/FileMaster.h>
#include <util/misc/ioUtil.h>

#include <sstream>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   MdWriteRestart::MdWriteRestart(MdSimulation& simulation) 
    : Diagnostic(),
      filename_(),
      simulationPtr_(&simulation)
   {  setClassName("MdWriteRestart"); }

   /*
   * Read interval and outputFileName. 
   */
   void MdWriteRestart::readParameters(std::istream& in) 
   {
      readInterval(in);
      read<std::string>(in, "fileName", filename_);
   }

   /*
   * Dump configuration to file
   */
   void MdWriteRestart::sample(long iStep) 
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
