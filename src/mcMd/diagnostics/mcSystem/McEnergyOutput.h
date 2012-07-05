#ifndef MCMD_MC_ENERGY_OUTPUT_H
#define MCMD_MC_ENERGY_OUTPUT_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/diagnostics/SystemDiagnostic.h> // base class template
#include <mcMd/mcSimulation/McSystem.h>        // base template parameter

namespace McMd
{

   using namespace Util;

   /**
   * Diagnostic to output total potential energy.
   *
   * \ingroup McMd_Diagnostic_Module
   */
   class McEnergyOutput : public SystemDiagnostic<McSystem>
   {

   public:
   
      /// Constructor.
      McEnergyOutput(McSystem& system);

      /// Read output file and nStepPerSample.
      virtual void readParam(std::istream& in);
 
      /// Evaluate energy and print.
      void sample(long iStep);

      /// Output final summary and file format
      virtual void output();

   private:

     /// Output file stream
     std::ofstream outputFile_;
   
   };

}
#endif 
