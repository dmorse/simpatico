#ifndef INTER_INTRA_LINK_H
#define INTER_INTRA_LINK_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/diagnostics/SystemDiagnostic.h>  // base class template
#include <mcMd/simulation/System.h>             // class template parameter
#include <util/accumulators/Average.h>          // member

#include <cstdio> 
#include <cstring> 

namespace McMd
{

   using namespace Util;

   class Species;

   /**
   * Number of inter and intramolecular links
   *
   *
   * This class calculates the average number of inter and intramolecular links by accumulating an
   * average value for the above quantities separately, and optionally outputs block
   * averages to file at an interval specified by the input parameter
   * nSamplePerBlock. No block averages are output if nSamplePerBlock = 0.
   *
   * \ingroup Diagnostic_Module
   */
   class InterIntraLink : public SystemDiagnostic<System>
   {
   
   public:

      /**
      * Constructor.
      *
      * \param system reference to parent System object
      */
      InterIntraLink(System &system);
   
      /**
      * Read parameters from file, and allocate data array.
      *
      * Input format:
      *
      *   - int    interval        : sampling interval
      *   - string outputFileName  : base name for output file(s)
      *   - int    nSamplePerBlock : interval for output of block averages
      *   - int    speciesId       : integer id for Species of interest
      *
      * No block averages are output if nSamplePerBlock = 0. Otherwise,
      * block averages are output to a file named (outputFileName).dat. 
      *
      * \param in parameter input stream
      */
      virtual void readParameters(std::istream& in);
   
      /** 
      * Clear accumulator.
      */
      virtual void setup();
   
      /**
      * Evaluate inter and intramolecule number of links, add to ensemble.
      *
      * \param iStep step counter
      */
      virtual void sample(long iStep);
   
      /**
      * Output results at end of simulation.
      */
      virtual void output();

   private:

      /// Output file stream
      std::ofstream outputFile_;

      /// Average object for intermolecular links  - statistical accumulator
      Average  accumulatorInter_;
      
      /// Average object for intramolecular links - statistical accumulator
      Average  accumulatorIntra_;  

      /// Number of samples per block average output.
      int nSamplePerBlock_;

      /// Pointer to relevant Species.
      int     speciesId_;
   
   };

}
#endif
