#ifndef MCMD_MD_ENERGY_OUTPUT_H
#define MCMD_MD_ENERGY_OUTPUT_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/diagnostics/SystemDiagnostic.h>  // base class template
#include <mcMd/mdSimulation/MdSystem.h>         // base template parameter
#include <util/global.h> 

namespace McMd
{

   using namespace Util;

   /**
   * Diagnostic to output total potential and kinetic energies.
   * 
   * \ingroup Diagnostic_Module
   */
   class MdEnergyOutput : public SystemDiagnostic<MdSystem>
   {

   public:
   
      /// Constructor.
      MdEnergyOutput(MdSystem& system);

      /// Read output file and nStepPerSample.
      virtual void readParam(std::istream& in);
 
      /// Evaluate energy and print.
      virtual void sample(long iStep);

      /// Output final summary and file format
      virtual void output();

   private:

      /// Output file stream
      std::ofstream outputFile_;
   
   };

}
#endif 
