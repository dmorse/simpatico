#ifndef DDMD_OUTPUT_ENERGY_CPP
#define DDMD_OUTPUT_ENERGY_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "OutputEnergy.h"
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
   OutputEnergy::OutputEnergy(System& system) 
    : Diagnostic(system),
      nSample_(0),
      isInitialized_(false)
   {}

   /*
   * Read interval and outputFileName. 
   */
   void OutputEnergy::readParam(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);
      isInitialized_ = true;
   }

   /*
   * Read interval and outputFileName. 
   */
   void OutputEnergy::setup() 
   {  
       nSample_ = 0; 
      //std::string filename;
      //filename  = outputFileName();
      //system().fileMaster().openOutputFile(filename, outputFile_);
   }

   /*
   * Dump configuration to file
   */
   void OutputEnergy::sample(long iStep) 
   {
      if (isAtInterval(iStep))  {
         System& sys = system();
         sys.computeKineticEnergy();
         sys.computePotentialEnergies();
         if (sys.domain().isMaster()) {
            double kinetic = sys.kineticEnergy();
            double potential = sys.potentialEnergy();
            std::cout << Int(iStep, 10)
                      << Dbl(kinetic, 20)
                      << Dbl(potential, 20)
                      << Dbl(kinetic + potential, 20)
                      << std::endl;
         }

         ++nSample_;
      }
   }

}
#endif 
