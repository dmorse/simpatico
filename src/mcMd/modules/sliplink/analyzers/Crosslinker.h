#ifndef CROSSLINKER_H
#define CROSSLINKER_H

/*
* MolMcD - Monte Carlo and Molecular Dynamics Simulator for Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/SystemAnalyzer.h>
#include <mcMd/simulation/System.h>
#include <mcMd/neighbor/CellList.h>

namespace McMd
{
    
   using namespace Util;
  
   /**
   * Analyzer to create crosslinks and output the resulting configuration.
   * 
   * \ingroup Analyzer_Module
   */
   class Crosslinker : public SystemAnalyzer<System>
   {

   public:
   
      /// Constructor.
      Crosslinker(System& system);

      /// Read output file and nStepPerSample.
      virtual void readParameters(std::istream& in);

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
