/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AverageAnalyzer.h"
#include <ddMd/simulation/Simulation.h>
#include <util/accumulators/Average.h>   
#include <util/format/Int.h>
#include <util/format/Dbl.h>
#include <util/misc/FileMaster.h>
#include <util/misc/ioUtil.h>

#include <sstream>

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   AverageAnalyzer::AverageAnalyzer(Simulation& simulation) 
    : Analyzer(simulation),
      AverageMixIn(simulation.fileMaster())
   {  setClassName("AverageAnalyzer"); }

   /*
   * Destructor.
   */
   AverageAnalyzer::~AverageAnalyzer() 
   {}

   /*
   * Read interval and outputFileName. 
   */
   void AverageAnalyzer::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);
      readNSamplePerBlock(in, *this);
      if (simulation().domain().isMaster()) {
         initializeAccumulator();
      }
   }

   /*
   * Load internal state from an archive.
   */
   void AverageAnalyzer::loadParameters(Serializable::IArchive &ar)
   {
      loadInterval(ar);
      loadOutputFileName(ar);
      loadNSamplePerBlock(ar, *this);
      if (simulation().domain().isMaster()) {
         loadAccumulator(ar);
      }
   }

   /*
   * Save internal state to an archive.
   */
   void AverageAnalyzer::save(Serializable::OArchive &ar)
   {
      UTIL_CHECK(simulation().domain().isMaster());
      UTIL_CHECK(hasAccumulator());

      saveInterval(ar);
      saveOutputFileName(ar);
      saveNSamplePerBlock(ar);
      saveAccumulator(ar);
   }

   /*
   * Clear accumulator (do nothing on slave processors).
   */
   void AverageAnalyzer::clear() 
   {   
      if (simulation().domain().isMaster()){ 
         clearAccumulator();
      }
   }
 
   /*
   * Open outputfile
   */ 
   void AverageAnalyzer::setup()
   {
      if (simulation().domain().isMaster()) {
         UTIL_CHECK(hasAccumulator());
         if (nSamplePerBlock() > 0) {
            openOutputFile(outputFileName(".dat"));
         }
      }
   }

   /*
   * Compute value.
   */
   void AverageAnalyzer::sample(long iStep) 
   {
      if (!isAtInterval(iStep)) {
         UTIL_THROW("Time step index is not a multiple of interval");
      }
      compute();
      if (simulation().domain().isMaster()) {
         updateAccumulator(iStep, interval());
      }
   }

   /*
   * Output results to several files after simulation is completed.
   */
   void AverageAnalyzer::output()
   {
      if (hasAccumulator()) {

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

}
