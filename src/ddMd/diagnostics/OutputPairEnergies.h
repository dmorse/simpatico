#ifndef DDMD_OUTPUT_PAIR_ENERGIES_H
#define DDMD_OUTPUT_PAIR_ENERGIES_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/diagnostics/Diagnostic.h>
#include <ddMd/simulation/Simulation.h>
#include <util/containers/DMatrix.h>

namespace DdMd
{

   using namespace Util;

   /**
   * Periodically write pair energies to file.
   *
   * \ingroup DdMd_Diagnostic_Module
   */
   class OutputPairEnergies : public Diagnostic
   {
   
   public:
   
      /**
      * Constructor.
      *
      * \param simulation parent Simulation object. 
      */
      OutputPairEnergies(Simulation& simulation);
   
      /**
      * Destructor.
      */
      virtual ~OutputPairEnergies()
      {} 
   
      /**
      * Read interval and output file name.
      *
      * \param in input parameter file
      */
      virtual void readParameters(std::istream& in);
   
      /**
      * Clear nSample counter.
      */
      virtual void clear();
  
      /**
      * Compute and output pair energies periodically.
      *
      * \param iStep MC step index
      */
      virtual void sample(long iStep);

   private:
 
      // Output file stream.
      std::ofstream outputFile_;

      /// Number of samples.
      long    nSample_;

      /// Has readParam been called?
      long    isInitialized_;
   
   };

}
#endif 
