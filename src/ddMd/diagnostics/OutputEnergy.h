#ifndef DDMD_OUTPUT_ENERGY_H
#define DDMD_OUTPUT_ENERGY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/diagnostics/Diagnostic.h>
#include <ddMd/system/System.h>

namespace DdMd
{

   using namespace Util;

   /**
   * Periodically write system energies to file.
   *
   * \ingroup Diagnostic_Module
   */
   class OutputEnergy : public Diagnostic
   {
   
   public:
   
      /**
      * Constructor.
      *
      * \param system parent System object. 
      */
      OutputEnergy(System& system);
   
      /**
      * Destructor.
      */
      virtual ~OutputEnergy()
      {} 
   
      /**
      * Read dumpPrefix and interval.
      *
      * \param in input parameter file
      */
      virtual void readParam(std::istream& in);
   
      /**
      * Clear nSample counter.
      */
      virtual void setup();
  
      /**
      * Dump configuration to file
      *
      * \param iStep MC step index
      */
      virtual void sample(long iStep);

   private:
 
      // Output file stream
      std::ofstream outputFile_;

      /// Number of configurations dumped thus far (first dump is zero).
      long    nSample_;
   
      /// Has readParam been called?
      long    isInitialized_;
   
   };

}
#endif 
