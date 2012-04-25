#ifndef SS_CHAIN_DIST_H
#define SS_CHAIN_DIST_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/diagnostics/SystemDiagnostic.h>    // base class template
#include <mcMd/simulation/System.h>               // base class template parameter
#include <util/accumulators/Distribution.h>
#include <util/containers/DArray.h>

#include <util/global.h>

namespace McMd
{

   using namespace Util;

   /**
   * SSChainDist evaluates the distribution of slip-springs along the chains.
   * Suggested values for the histogram:
   * min   -1.0
   * max    chainlength                        
   * nBin   chainlength+1
   *
   * \ingroup Diagnostic_Module
   */
   class SSChainDist : public SystemDiagnostic<System>
   {
   
   public:

      /**	
      * Constructor.
      *
      * \param system reference to parent System object
      */
      SSChainDist(System &system);

      /**
      * Read parameters from file.  
      *
      * \param in input parameter stream
      */
      virtual void readParam(std::istream& in);

      /** 
      * Clear accumulator.
      */
      virtual void setup();
   
      /** 
      * Add atoms attached to links to SSChainDist histogram.
      *
      * \param iStep step counter
      */
      void sample(long iStep);

      /** 
      * Output results to output file.
      */
      virtual void output();

   private:

      // Output file stream
      std::ofstream outputFile_;

      // Distribution statistical accumulator
      Distribution  accumulator_;

   };

}
#endif
