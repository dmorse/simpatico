#ifndef DDMD_OUTPUT_PRESSURE_CPP
#define DDMD_OUTPUT_PRESSURE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "OutputPressure.h"
//#include <ddMd/util/FileMaster.h>
#include <util/util/ioUtil.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

#include <sstream>

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   OutputPressure::OutputPressure(Simulation& simulation) 
    : Diagnostic(simulation),
      nSample_(0),
      isInitialized_(false)
   {}

   /*
   * Read interval and outputFileName. 
   */
   void OutputPressure::readParam(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);
      isInitialized_ = true;
   }

   /*
   * Read interval and outputFileName. 
   */
   void OutputPressure::setup() 
   {  
      nSample_ = 0; 
      std::string filename;
      filename  = outputFileName();
      simulation().fileMaster().openOutputFile(filename, outputFile_);
   }

   /*
   * Dump configuration to file
   */
   void OutputPressure::sample(long iStep) 
   {
      if (isAtInterval(iStep))  {
         Simulation& sys = simulation();
         sys.computeVirialStress();
         sys.computeKineticStress();
         if (sys.domain().isMaster()) {
            double virial  = sys.virialPressure();
            double kinetic = sys.kineticPressure();
            outputFile_ << Int(iStep, 10)
                        << Dbl(kinetic, 20)
                        << Dbl(virial, 20)
                        << Dbl(kinetic + virial, 20)
                        << std::endl;
         }

         ++nSample_;
      }
   }

}
#endif 
