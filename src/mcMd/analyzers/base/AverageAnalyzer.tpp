#ifndef MCMD_AVERAGE_ANALYZER_TPP
#define MCMD_AVERAGE_ANALYZER_TPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AverageAnalyzer.h"
#include <util/accumulators/Average.h>

namespace McMd
{

   using namespace Util;
   using namespace Simp;

   /*
   * Constructor.
   */
   template <class SystemType>
   AverageAnalyzer<SystemType>::AverageAnalyzer(SystemType& system) 
    : SystemAnalyzer<SystemType>(system),
      AverageMixIn(system.fileMaster())
   {  setClassName("AverageAnalyzer"); }

   /*
   * Destructor.
   */
   template <class SystemType>
   AverageAnalyzer<SystemType>::~AverageAnalyzer() 
   {}

   /*
   * Read interval and outputFileName. 
   */
   template <class SystemType>
   void AverageAnalyzer<SystemType>::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);
      readNSamplePerBlock(in, *this);
      initializeAccumulator();
      UTIL_CHECK(nSamplePerBlock() >= 0);
      UTIL_CHECK(accumulator().nSamplePerBlock() == nSamplePerBlock());
      if (nSamplePerBlock() > 0) {
         openOutputFile(outputFileName(".dat"));
      }
   }

   /*
   * Load internal state from an archive.
   */
   template <class SystemType>
   void AverageAnalyzer<SystemType>::loadParameters(Serializable::IArchive &ar)
   {
      loadInterval(ar);
      loadOutputFileName(ar);
      loadNSamplePerBlock(ar, *this);
      loadAccumulator(ar);
      UTIL_CHECK(nSamplePerBlock() >= 0);
      UTIL_CHECK(accumulator().nSamplePerBlock() == nSamplePerBlock());
      if (nSamplePerBlock() > 0) {
         openOutputFile(outputFileName(".dat"));
      }
   }

   /*
   * Save internal state to an archive.
   */
   template <class SystemType>
   void AverageAnalyzer<SystemType>::save(Serializable::OArchive &ar)
   {
      UTIL_CHECK(hasAccumulator());

      SystemAnalyzer<SystemType>::save(ar);
      saveNSamplePerBlock(ar);
      saveAccumulator(ar);
   }

   /*
   * Clear accumulator.
   */
   template <class SystemType>
   void AverageAnalyzer<SystemType>::clear() 
   {   
      clearAccumulator();
   }
 
   /*
   * Open outputfile
   */ 
   template <class SystemType>
   void AverageAnalyzer<SystemType>::setup()
   {
      UTIL_CHECK(hasAccumulator());
      clearAccumulator();
   }

   /*
   * Get current value (default implementation).
   */ 
   template <class SystemType>
   double AverageAnalyzer<SystemType>::value()
   {  return value_; }

   /*
   * Compute current value and update accumulator and output.
   */
   template <class SystemType>
   void AverageAnalyzer<SystemType>::sample(long iStep) 
   {
      if (isAtInterval(iStep)) {
         compute();
         updateAccumulator(iStep, interval());
      }
   }

   /*
   * Output results to several files after simulation is completed.
   */
   template <class SystemType>
   void AverageAnalyzer<SystemType>::output()
   {
      UTIL_CHECK(hasAccumulator());

      // Close data (*.dat) file, if any
      if (outputFile().is_open()) {
         outputFile().close();
      }

      // Write parameter (*.prm) file
      openOutputFile(outputFileName(".prm"));
      ParamComposite::writeParam(outputFile());
      outputFile().close();

      // Write average (*.ave) and error analysis (*.aer) files
      outputAccumulator(outputFileName());
   }

}
#endif
