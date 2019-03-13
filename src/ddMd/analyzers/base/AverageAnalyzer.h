#ifndef DDMD_AVERAGE_ANALYZER_H
#define DDMD_AVERAGE_ANALYZER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/analyzers/Analyzer.h>

namespace Util { 
   class Average;
   class FileMaster;
}

namespace DdMd
{

   class Simulation;

   using namespace Util;

   /**
   * Analyze average and block averages of a single floating point variable.
   *
   * This class evaluates the average of a sampled float point variable, and
   * optionally writes block averages to a data file during a simulation. 
   * It should be used as a base class of Analyzers that evaluate averages
   * and (optionally) output block averages for specific physical variables.
   *
   * \ingroup DdMd_Analyzer_Base_Module
   */
   class AverageAnalyzer : public Analyzer
   {
   
   public:
   
      /**
      * Constructor.
      *
      * \param simulation  parent Simulation object. 
      */
      AverageAnalyzer(Simulation& simulation);
   
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

      /**
      * Does this processor have an Average accumulator?
      */
      bool hasAccumulator() const;

      /**
      * Get Average accumulator.
      *
      * Call only on master.
      *
      * \param i integer index of value.
      */
      const Average& accumulator() const;

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

      /**
      * Access output file by reference.
      */
      std::ofstream& outputFile();

      /**
      * Open an output file. 
      */
      void openOutputFile(std::string filename, std::ofstream& file);

   private:

      /// Output file stream.
      std::ofstream  outputFile_;

      /// Pointer to Average object (only instantiated on master processor)
      Average *accumulatorPtr_;

      /// Pointer to a FileMaster
      FileMaster* fileMasterPtr_;

      /// Number of samples per block average output.
      int nSamplePerBlock_;
   
      /// Has readParam been called?
      bool isInitialized_;
   
   };


   // Inline functions

   /*
   * Does this processor have an accumulator?
   */
   inline
   bool AverageAnalyzer::hasAccumulator() const
   {  return (bool)(accumulatorPtr_); }

   /*
   * Get accumulator associated with a variable.
   */
   inline
   const Average& AverageAnalyzer::accumulator() const
   {
      UTIL_CHECK(accumulatorPtr_);
      return *accumulatorPtr_;
   }

   /*
   * Access output file by reference.
   */
   inline
   std::ofstream& AverageAnalyzer::outputFile()
   {  return outputFile_; }

}
#endif 
