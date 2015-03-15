#ifndef DDMD_AVERAGE_ANALYZER_H
#define DDMD_AVERAGE_ANALYZER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/analyzers/Analyzer.h>

namespace Util { 
   class Average;
}

namespace DdMd
{

   class Simulation;

   using namespace Util;

   /**
   * Base class for analyzers that evaluate averages of global scalars.
   *
   * This class evaluates the average of a sampled float point values,
   * and optionally writes block averages to a data file during the run.
   *
   * \ingroup DdMd_Analyzer_Module
   */
   class AverageAnalyzer : public Analyzer
   {
   
   public:
   
      /**
      * Constructor.
      *
      * \param simulation parent Simulation object. 
      */
      AverageAnalyzer(Simulation& simulation);
   
      /**
      * Destructor.
      */
      virtual ~AverageAnalyzer(); 
   
      /**
      * Read dumpPrefix and interval.
      *
      * \param in input parameter file
      */
      virtual void readParameters(std::istream& in);
   
      /**
      * Load internal state from an archive.
      *
      * \param ar input/loading archive
      */
      virtual void loadParameters(Serializable::IArchive &ar);

      /**
      * Save internal state to an archive.
      *
      * \param ar output/saving archive
      */
      virtual void save(Serializable::OArchive &ar);
  
      /**
      * Clear accumulator on master, do nothing on other processors.
      */
      virtual void clear();
  
      /**
      * Setup before loop - open output file.
      */
      virtual void setup();

      /**
      * Compute a sampled value and add it to a sequence.
      *
      * \param iStep MD step index
      */
      virtual void sample(long iStep);

      /**
      * Dump configuration to file
      */
      virtual void output();

   protected:

      /**
      * Compute value of sampled quantity.
      *
      * Call on all processors.
      */
      virtual void compute() = 0;

      /**
      * Get current value, set by compute function.
      *
      * Call only on master.
      */
      virtual double value() = 0;

   private:

      /// Output file stream.
      std::ofstream  outputFile_;

      /// Pointer to Average object (only instantiated on master processor)
      Average *accumulatorPtr_;

      /// Number of samples per block average output.
      int nSamplePerBlock_;
   
      /// Has readParam been called?
      bool isInitialized_;
   
   };

}
#endif 
