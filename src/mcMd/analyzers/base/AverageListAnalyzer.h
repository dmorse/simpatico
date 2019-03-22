#ifndef MCMD_AVERAGE_LIST_ANALYZER_H
#define MCMD_AVERAGE_LIST_ANALYZER_H

/*
* Simpatico - SystemType Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/base/SystemAnalyzer.h>
#include <simp/analysis/AverageListMixIn.h>

namespace McMd
{

   using namespace Util;
   using namespace Simp;

   /**
   * Analyze averages and block averages of several real variables.
   *
   * This class evaluates the average of several sampled real variables, and
   * optionally writes block averages to a data file during a simulation. 
   * It is intended for use as a base class for Analyzers that evaluate 
   * averages and (optionally) block averages for specific physical variables.
   *
   * \ingroup McMd_Analyzer_Base_Module
   */
   template <class SystemType>
   class AverageListAnalyzer : public SystemAnalyzer<SystemType>, 
                               public AverageListMixIn
   {
   
   public:
   
      /**
      * Constructor.
      *
      * \param system parent SystemType object. 
      */
      AverageListAnalyzer(SystemType& system);
   
      /**
      * Destructor.
      */
      virtual ~AverageListAnalyzer(); 

      /**
      * Read interval, outputFileName and (optionally) nSamplePerBlock.
      *
      * The optional variable nSamplePerBlock defaults to 0, which disables
      * computation and output of block averages. Setting nSamplePerBlock = 1
      * outputs every sampled value. 
      *
      * \param in  input parameter file
      */
      virtual void readParameters(std::istream& in);
   
      /**
      * Load internal state from an input archive.
      *
      * \param ar  input/loading archive
      */
      virtual void loadParameters(Serializable::IArchive &ar);

      /**
      * Save internal state to an output archive.
      *
      * \param ar  output/saving archive
      */
      virtual void save(Serializable::OArchive &ar);
  
      /**
      * Clear accumulator on master, do nothing on other processors.
      */
      virtual void clear();
  
      /**
      * Setup before loop. Opens an output file, if any.
      */
      virtual void setup();

      /**
      * Compute a sampled value and update the accumulator.
      *
      * \param iStep  MD time step index
      */
      virtual void sample(long iStep);

      /**
      * Write final results to file after a simulation.
      */
      virtual void output();
   
      using SystemAnalyzer<SystemType>::interval;
      using SystemAnalyzer<SystemType>::isAtInterval;
      using SystemAnalyzer<SystemType>::outputFileName;

   protected:

      using SystemAnalyzer<SystemType>::setClassName;
      using SystemAnalyzer<SystemType>::readInterval;
      using SystemAnalyzer<SystemType>::readOutputFileName;
      using SystemAnalyzer<SystemType>::loadInterval;
      using SystemAnalyzer<SystemType>::loadOutputFileName;
      using SystemAnalyzer<SystemType>::save;

      using AverageListMixIn::outputFile;
      using AverageListMixIn::openOutputFile;

      using AverageListMixIn::readNSamplePerBlock;
      using AverageListMixIn::loadNSamplePerBlock;
      using AverageListMixIn::loadAccumulators;
      using AverageListMixIn::saveNSamplePerBlock;
      using AverageListMixIn::saveAccumulators;

   };

}
#endif 
