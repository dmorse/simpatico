#ifndef MCMD_AVERAGE_ANALYZER_H
#define MCMD_AVERAGE_ANALYZER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/base/SystemAnalyzer.h>    // base class
#include <simp/analysis/AverageMixIn.h>            // base class

namespace McMd
{

   using namespace Util;
   using namespace Simp;

   /**
   * Analyze average and block averages of a single floating point variable.
   *
   * This class evaluates the average of a sampled float point variable, and
   * optionally writes block averages to a data file during a simulation. 
   * It should be used as a base class of Analyzers that evaluate averages
   * and (optionally) output block averages for specific physical variables.
   *
   * \ingroup McMd_Analyzer_Base_Module
   */
   template <class SystemType>
   class AverageAnalyzer : public SystemAnalyzer<SystemType>, 
                           public AverageMixIn
   {

   public:
   
      /**
      * Constructor.
      *
      * \param system parent system object. 
      */
      AverageAnalyzer(SystemType& system);
   
      /**
      * Destructor.
      */
      virtual ~AverageAnalyzer(); 
   
      /**
      * Read interval, outputFileName and (optionally) nSamplePerBlock.
      *
      * The optional variable nSamplePerBlock defaults to 0, which disables
      * computation and output of block averages. Set nSamplePerBlock = 1
      * to output every sampled value. 
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

      // Current value, should be set by compute function.
      double value_;
  
      /**
      * Get current value (default implementation)
      */ 
      virtual double value();

      using SystemAnalyzer<SystemType>::setClassName;
      using SystemAnalyzer<SystemType>::readInterval;
      using SystemAnalyzer<SystemType>::readOutputFileName;
      using SystemAnalyzer<SystemType>::loadInterval;
      using SystemAnalyzer<SystemType>::loadOutputFileName;
      using SystemAnalyzer<SystemType>::save;

      using AverageMixIn::outputFile;
      using AverageMixIn::openOutputFile;

      using AverageMixIn::readNSamplePerBlock;
      using AverageMixIn::loadNSamplePerBlock;
      using AverageMixIn::loadAccumulator;
      using AverageMixIn::saveNSamplePerBlock;
      using AverageMixIn::saveAccumulator;

   };

}
#endif 
