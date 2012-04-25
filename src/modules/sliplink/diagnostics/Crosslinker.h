#ifndef CROSSLINKER_H
#define CROSSLINKER_H

/*
* MolMcD - Monte Carlo and Molecular Dynamics Simulator for Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/diagnostics/SystemDiagnostic.h>
#include <mcMd/simulation/System.h>
#include <mcMd/neighbor/CellList.h>

namespace McMd
{
    
   using namespace Util;
  
   /**
   * Diagnostic to create crosslinks and output the resulting configuration.
   * 
   * \ingroup Diagnostic_Module
   */
   class Crosslinker : public SystemDiagnostic<System>
   {

   public:
   
      /// Constructor.
      Crosslinker(System& system);

      /// Read output file and nStepPerSample.
      virtual void readParam(std::istream& in);

      void setup();
 
      /// Create crosslinks and dump configuration to file.
      void sample(long iStep);

   private:

      /// Output file stream
      std::ofstream outputFile_;

      /// Number of configurations dumped thus far (first dump is zero).
      long    nSample_;

      /// Private CellList, used to search for link neighbors.
      CellList cellList_;

      double cutoff_;
      double probability_;
   
   };

}

#endif
