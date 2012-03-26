#ifndef MCMD_SEMI_GRAND_DISTRIBUTION_H
#define MCMD_SEMI_GRAND_DISTRIBUTION_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/diagnostics/SystemDiagnostic.h>  // base class template
#include <mcMd/mcSimulation/McSystem.h>         // base template parameter
#include <util/accumulators/IntDistribution.h>  // member

namespace McMd
{

   using namespace Util;
   class Species;
   class SpeciesMutator;

   /**
   * Calculate distribution of type indices for mutable species.
   *
   * \ingroup McMd_Diagnostic_Module
   */
   class SemiGrandDistribution : public SystemDiagnostic<McSystem>
   {

   public:
   
      /// Constructor.
      SemiGrandDistribution(McSystem& system);

      /// Read output file and nStepPerSample.
      virtual void readParam(std::istream& in);
 
      /// Evaluate energy and print.
      void sample(long iStep);

      /// Output final summary and file format
      virtual void output();

   private:

     /// Output file stream.
     std::ofstream outputFile_;

     /// Histogram of values.
     IntDistribution distribution_;

     /// Species index.
     int speciesId_;

     /// Maximum possible number of molecules. 
     int moleculeCapacity_;

     /// Pointer to Species.
     Species* speciesPtr_;

     /// Pointer to associated SpeciesMutator.
     SpeciesMutator* mutatorPtr_;
   
   };

}
#endif 
