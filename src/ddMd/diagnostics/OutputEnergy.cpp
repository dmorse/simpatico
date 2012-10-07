#ifndef DDMD_OUTPUT_ENERGY_CPP
#define DDMD_OUTPUT_ENERGY_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "OutputEnergy.h"
//#include <util/misc/FileMaster.h>
#include <util/misc/ioUtil.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

#include <sstream>

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   OutputEnergy::OutputEnergy(Simulation& simulation) 
    : Diagnostic(simulation),
      nSample_(0),
      isInitialized_(false)
   {}

   /*
   * Read interval and outputFileName. 
   */
   void OutputEnergy::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);
      isInitialized_ = true;
   }

   /*
   * Read interval and outputFileName. 
   */
   void OutputEnergy::clear() 
   {  nSample_ = 0;  }

   /*
   * Dump configuration to file
   */
   void OutputEnergy::sample(long iStep) 
   {
      if (isAtInterval(iStep))  {
         Simulation& sys = simulation();
         sys.computeKineticEnergy();
         sys.computePotentialEnergies();
         if (sys.domain().isMaster()) {
            double kinetic = sys.kineticEnergy();
            double potential = sys.potentialEnergy();
            Log::file() << Int(iStep, 10)
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
