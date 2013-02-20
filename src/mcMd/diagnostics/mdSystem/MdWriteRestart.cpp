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
//#include <util/misc/FileMaster.h>
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
      simulationPtr_(&simulation)
   {  setClassName("MdWriteRestart"); }

   /*
   * Read interval and outputFileName. 
   */
   void MdWriteRestart::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);
   }

   // Using default Diagnostic::save and Diagnostic::load.

   /*
   * Write a restart file. 
   */
   void MdWriteRestart::sample(long iStep) 
   {
      if (isAtInterval(iStep))  {
         simulationPtr_->writeRestart(outputFileName());
      }
   }
  
}
#endif 
