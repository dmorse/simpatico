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
   
      // Analyzer interface.

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

      // Extension of public interface.

      /**
      * Get value of nSamplePerBlock.
      * 
      * If nSamplePerBlock == 0, output of block averages is disabled.
      * For nSamplePerBlock > 0, the value is the number of sampled values
      * averaged in each block. 
      */
      int nSamplePerBlock() const;

      /**
      * Does this processor have an Average accumulator?
      */
      bool hasAccumulator() const;

      /**
      * Get Average accumulator.
      *
      * \pre hasAccumulator() == true
      *
      * \param i integer index of value.
      */
      const Average& accumulator() const;

   protected:

      /**
      * Read nSamplePerBlock parameter from file.
      *
      * \param in parameter file.
      */ 
      void readNSamplePerBlock(std::istream& in);

      /**
      * Instantiate a new Average accumulator and set nSamplePerBlock.
      *
      * \pre hasAccumulator()
      * \pre nSamplePerBlock is set and positive
      */ 
      void initializeAccumulator();

      /**
      * Clear internal state of the accumulator.
      *
      * \pre hasAccumulator()
      */ 
      void clearAccumulator();

      /**
      * Load nSamplePerBlock parameter from an archive.
      *
      * \param ar input archive
      */ 
      void loadNSamplePerBlock(Serializable::IArchive &ar);

      /**
      * Instantiate an accumulator and load data from an archive.
      *
      * \param ar input archive
      */ 
      void loadAccumulator(Serializable::IArchive &ar);

      /**
      * Save nSamplePerBlock parameter to an archive.
      *
      * \param ar output archive
      */ 
      void saveNSamplePerBlock(Serializable::OArchive &ar);

      /**
      * Save accumulator to an archive.
      *
      * \param ar output archive
      */ 
      void saveAccumulator(Serializable::OArchive &ar);

      /**
      * Compute value of sampled quantity.
      *
      * Call on all processors.
      */
      virtual void compute() = 0;

      /**
      * Get current value, set by compute function.
      *
      * \pre hasAccumulator() == true
      */
      virtual double value() = 0;

      /**
      * Add current value to accumulator, output block average if needed.
      *
      * \pre hasAccumulator() == true
      * \param iStep simulation step counter
      */
      void updateAccumulator(long iStep);

      /**
      * Write results of statistical analysis to files.
      *
      * \pre hasAccumulator() == true
      */
      void outputAccumulator();

      /**
      * Access output file by reference.
      */
      std::ofstream& outputFile();

      /**
      * Open the output file. 
      *
      * \param filename base file name, without output prefix
      */
      void openOutputFile(std::string filename);

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
   int AverageAnalyzer::nSamplePerBlock() const
   {  return nSamplePerBlock_; }

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
