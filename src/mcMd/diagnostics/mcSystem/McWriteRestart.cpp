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
//#include <util/misc/FileMaster.h>
#include <util/misc/ioUtil.h>

#include <sstream>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   McWriteRestart::McWriteRestart(McSimulation& simulation) 
    : Diagnostic(),
      simulationPtr_(&simulation)
   {  setClassName("McWriteRestart"); }

   /*
   * Read interval and outputFileName. 
   */
   void McWriteRestart::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);
   }

   // Using Diagnostic::loadParameters and Diagnostic::save

   /*
   * Write state to restart file.
   */
   void McWriteRestart::sample(long iStep) 
   {
      if (isAtInterval(iStep))  {
         simulationPtr_->writeRestart(outputFileName());
      }
   }
  
}
#endif 
