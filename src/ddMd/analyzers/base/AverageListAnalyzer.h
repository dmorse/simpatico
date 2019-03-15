#ifndef DDMD_AVERAGE_LIST_ANALYZER_H
#define DDMD_AVERAGE_LIST_ANALYZER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/analyzers/Analyzer.h>
#include <util/containers/DArray.h>
#include <string>

namespace Util {
   class Average;
   class FileMaster;
}

namespace DdMd
{

   class Simulation;

   using namespace Util;

   /**
   * Analyze averages and block averages of several real variables.
   *
   * This class evaluates the average of several sampled real variables, and
   * optionally writes block averages to a data file during a simulation. 
   * It is intended for use as a base class for Analyzers that evaluate 
   * averages and (optionally) block averages for specific physical variables.
   *
   * \ingroup DdMd_Analyzer_Base_Module
   */
   class AverageListAnalyzer : public Analyzer
   {
   
   public:
   
      /**
      * Constructor.
      *
      * \param simulation  parent Simulation object. 
      */
      AverageListAnalyzer(Simulation& simulation);
   
      /**
      * Destructor.
      */
      virtual ~AverageListAnalyzer(); 

      // Member functions declared in Analyzer
   
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

      // Specialized public interface (not required by Analyzer base)

      /**
      * Get value of nSamplePer Block.
      */
      int nSamplePerBlock() const;

      /**
      * Get number of variables.
      *
      * Call only on processors that have accumulators.
      */
      int nValue() const;

      /**
      * Get name associated with value.
      *
      * Call only on processors that have accumulators.
      *
      * \param i integer index of name/value pair.
      */
      const std::string& name(int i) const;

      /**
      * Does this processor have accumulators?
      */
      bool hasAccumulators() const;

      /**
      * Get Average accumulator for a specific value.
      *
      * Call only on processors that have accumulators.
      *
      * \param i integer index of value.
      */
      const Average& accumulator(int i) const;

   protected:

      /**
      * Allocate accumulators and names arrays and set nValue.
      */
      void initializeAccumulators(int nValue);

      /**
      * Clear all accumulators.
      */
      void clearAccumulators();

      /**
      * Set name of variable. 
      *
      * Call only on master.
      *
      * \param i integer index of variable
      * \param name name of variable number i
      */
      void setName(int i, std::string name);

      /**
      * Read nSamplePerBlock parameter from file.
      *
      * \param in parameter file.
      */ 
      void readNSamplePerBlock(std::istream& in);

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
      void loadAccumulators(Serializable::IArchive &ar);

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
      void saveAccumulators(Serializable::OArchive &ar);

      /**
      * Set current value, used by compute function.
      *
      * \param i integer index of variable
      * \param value current value of variable
      */
      void setValue(int i, double value);

      /**
      * Compute current values of all variables.
      *
      * Call on all processors.
      */
      virtual void compute() = 0;

      /**
      * Get current value of a specific variable.
      *
      * Call only on master.
      *
      * \param i integer index of variable.
      */
      double value(int i) const;

      /**
      * Add current values to accumulators, output any block averages.
      *
      * \param iStep simulation step counter
      */
      void updateAccumulators(long iStep);

      /**
      * Write final accumulator data to files .
      */
      void outputAccumulators();

      /**
      * Access output file by references.
      */
      std::ofstream& outputFile();

      /**
      * Open an output file. 
      */
      void openOutputFile(std::string filename);

   private:

      /// Output file stream.
      std::ofstream  outputFile_;

      /// Array of Average objects (only allocated on master processor)
      DArray<Average> accumulators_;

      /// Array of current values (only allocated on master processor)
      DArray<double> values_;

      /// Array of value names (only allocated on master processor)
      DArray<std::string> names_;

      /// Pointer to a FileMaster
      FileMaster* fileMasterPtr_;

      /// Number of samples per block average output.
      int nSamplePerBlock_;

      /// Number of values.
      int nValue_;
 
      /// Does this processor have accumulators ?
      bool hasAccumulators_;
   
      /// Has readParam been called?
      bool isInitialized_;
   
   };

   // Inline functions

   /*
   * Get nSamplePerBlock.
   */
   inline
   int AverageListAnalyzer::nSamplePerBlock() const
   {  
      return nSamplePerBlock_; 
   }

   /*
   * Does this processor have accumulators?
   */
   inline
   bool AverageListAnalyzer::hasAccumulators() const
   {  return hasAccumulators_; }

   /*
   * Get nValue (number of variables).
   */
   inline
   int AverageListAnalyzer::nValue() const
   {  
      UTIL_CHECK(hasAccumulators());
      return nValue_; 
   }

   /*
   * Get current value of a variable, set by compute function.
   */
   inline
   double AverageListAnalyzer::value(int i) const
   {
      UTIL_CHECK(hasAccumulators());
      UTIL_CHECK(i >= 0 && i < nValue_);
      return values_[i];
   }

   /*
   * Get name of specific variable.
   */
   inline
   const std::string& AverageListAnalyzer::name(int i) const
   {
      UTIL_CHECK(hasAccumulators());
      UTIL_CHECK(i >= 0 && i < nValue_);
      return names_[i];
   }

   /*
   * Get accumulator associated with a variable.
   */
   inline
   const Average& AverageListAnalyzer::accumulator(int i) const
   {
      UTIL_CHECK(hasAccumulators());
      UTIL_CHECK(i >= 0 && i < nValue_);
      return accumulators_[i];
   }

   /*
   * Set current value of a variable.
   */
   inline
   void AverageListAnalyzer::setValue(int i, double value)
   {
      UTIL_CHECK(hasAccumulators());
      UTIL_CHECK(i >= 0 && i < nValue_);
      values_[i] = value;
   }

   /*
   * Access output file by reference.
   */
   inline
   std::ofstream& AverageListAnalyzer::outputFile()
   {  return outputFile_; }

}
#endif 
