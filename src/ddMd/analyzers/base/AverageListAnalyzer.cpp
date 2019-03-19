/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AverageListAnalyzer.h"
#include <ddMd/simulation/Simulation.h>

namespace DdMd
{

   using namespace Util;
   using namespace Simp;

   /*
   * Constructor.
   */
   AverageListAnalyzer::AverageListAnalyzer(Simulation& simulation) 
    : Analyzer(simulation),
      AverageListMixIn(simulation.fileMaster())
   {  setClassName("AverageListAnalyzer"); }

   /*
   * Destructor.
   */
   AverageListAnalyzer::~AverageListAnalyzer() 
   {}

   /*
   * Read interval and outputFileName. 
   */
   void AverageListAnalyzer::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);
      readNSamplePerBlock(in, *this);

      // Note: ReadParameters method of derived classes should call this,
      // determine nValue and then call initializeAccumulators(nValue).
   }

   /*
   * Load internal state from an archive.
   */
   void AverageListAnalyzer::loadParameters(Serializable::IArchive &ar)
   {
      loadInterval(ar);
      loadOutputFileName(ar);
      loadNSamplePerBlock(ar, *this);
      if (simulation().domain().isMaster()) {
         loadAccumulators(ar);
      }
   }

   /*
   * Save internal state to an archive.
   */
   void AverageListAnalyzer::save(Serializable::OArchive &ar)
   {
      UTIL_CHECK(simulation().domain().isMaster());
      UTIL_CHECK(hasAccumulators());
      saveInterval(ar);
      saveOutputFileName(ar);
      saveNSamplePerBlock(ar);
      saveAccumulators(ar);
   }

   /*
   * Clear accumulators (do nothing on slave processors).
   */
   void AverageListAnalyzer::clear() 
   {
      if (hasAccumulators()) {
         clearAccumulators();
      }
   }

   /*
   * Setup before simulation.
   */ 
   void AverageListAnalyzer::setup()
   {
      if (hasAccumulators()) {
         if (nSamplePerBlock()) {
            openOutputFile(outputFileName(".dat"));
         }
      }
   }

   /*
   * Compute and sample current values.
   */
   void AverageListAnalyzer::sample(long iStep) 
   {
      if (!isAtInterval(iStep)) {
         UTIL_THROW("Time step index is not a multiple of interval");
      }
      compute();
      if (hasAccumulators()) {
         updateAccumulators(iStep, interval());
      }
   }

   /*
   * Output results after a simulation is completed.
   */
   void AverageListAnalyzer::output()
   {
      if (hasAccumulators()) {

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

}
