#ifndef MCMD_AVERAGE_LIST_ANALYZER_TPP
#define MCMD_AVERAGE_LIST_ANALYZER_TPP

/*
* Simpatico - SystemType Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AverageListAnalyzer.h"

namespace McMd
{

   using namespace Util;
   using namespace Simp;

   /*
   * Constructor.
   */
   template <class SystemType>
   AverageListAnalyzer<SystemType>::AverageListAnalyzer(SystemType& system) 
    : SystemAnalyzer<SystemType>(system),
      AverageListMixIn(system.fileMaster())
   {  setClassName("AverageListAnalyzer"); }

   /*
   * Destructor.
   */
   template <class SystemType>
   AverageListAnalyzer<SystemType>::~AverageListAnalyzer() 
   {}

   /*
   * Read interval and outputFileName. 
   */
   template <class SystemType>
   void AverageListAnalyzer<SystemType>::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);
      readNSamplePerBlock(in, *this);
      if (nSamplePerBlock()) {
         openOutputFile(outputFileName(".dat"));
      }

      // Note: ReadParameters method of derived classes should call this,
      // determine nValue and then call initializeAccumulators(nValue).
   }

   /*
   * Load internal state from an archive.
   */
   template <class SystemType>
   void 
   AverageListAnalyzer<SystemType>::loadParameters(Serializable::IArchive &ar)
   {
      loadInterval(ar);
      loadOutputFileName(ar);
      loadNSamplePerBlock(ar, *this);
      loadAccumulators(ar);
      if (nSamplePerBlock()) {
         openOutputFile(outputFileName(".dat"));
      }
   }

   /*
   * Save internal state to an archive.
   */
   template <class SystemType>
   void AverageListAnalyzer<SystemType>::save(Serializable::OArchive &ar)
   {
      UTIL_CHECK(hasAccumulators());
      SystemAnalyzer<SystemType>::save(ar);
      saveNSamplePerBlock(ar);
      saveAccumulators(ar);
   }

   /*
   * Clear accumulators (do nothing on slave processors).
   */
   template <class SystemType>
   void AverageListAnalyzer<SystemType>::clear() 
   {
      UTIL_CHECK(hasAccumulators());
      clearAccumulators();
   }

   /*
   * Setup before system.
   */ 
   template <class SystemType>
   void AverageListAnalyzer<SystemType>::setup()
   {
      UTIL_CHECK(hasAccumulators());
      clear();
   }

   /*
   * Compute and sample current values.
   */
   template <class SystemType>
   void AverageListAnalyzer<SystemType>::sample(long iStep) 
   {
      UTIL_CHECK(hasAccumulators());
      if (!isAtInterval(iStep)) return;
      compute();
      updateAccumulators(iStep, interval());
   }

   /*
   * Output results after a system is completed.
   */
   template <class SystemType>
   void AverageListAnalyzer<SystemType>::output()
   {
      UTIL_CHECK(hasAccumulators());

      // Close data (*.dat) file, if any
      if (outputFile().is_open()) {
         outputFile().close();
      }

      // Write parameter (*.prm) file
      openOutputFile(outputFileName(".prm"));
      ParamComposite::writeParam(outputFile());
      outputFile().close();

      // Write average (*.ave) and error analysis (*.aer) files
      outputAccumulators(outputFileName());

   }

}
#endif
