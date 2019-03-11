#ifndef DDMD_AVERAGES_ANALYZER_H
#define DDMD_AVERAGES_ANALYZER_H

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
}

namespace DdMd
{

   class Simulation;

   using namespace Util;

   /**
   * Analyze multiple averages and block averages of real variables.
   *
   * This class evaluates the average of several sampled float point variable, 
   * and optionally writes block averages to a data file during a simulation. 
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

      /**
      * Get current values, set by compute function.
      *
      * Call only on master.
      *
      * \param i integer index of value.
      */
      double& value(int i);

      /**
      * Get name associated with value.
      *
      * Call only on master.
      *
      * \param i integer index of name/value pair.
      */
      const std::string& name(int i) const;

      /**
      * Get number of values.
      *
      * Call only on master.
      */
      int nValue() const;

   protected:

      /**
      * Set number of values.
      */
      void setNValue(int nValue);

      /**
      * Compute values of sampled quantities.
      *
      * Call on all processors.
      */
      virtual void compute() = 0;

   private:

      /// Output file stream.
      std::ofstream  outputFile_;

      /// Array of value names (only allocated on master processor)
      DArray<std::string> names_;

      /// Array of Average objects (only allocated on master processor)
      DArray<Average> accumulators_;

      /// Array of current values (only allocated on master processor)
      DArray<double> values_;

      /// Number of samples per block average output.
      int nSamplePerBlock_;

      /// Number of values.  
      int nValue_;
 
      /// Has readParam been called?
      bool isInitialized_;
   
   };

   // Inline functions

   /*
   * Get current value, set by compute function.
   */
   inline
   double& AverageListAnalyzer::value(int i)
   {
      UTIL_CHECK(nValue_ > 0);
      return values_[i];
   }

   /*
   * Get name of field.
   */
   inline
   const std::string& AverageListAnalyzer::name(int i) const
   {
      UTIL_CHECK(nValue_ > 0);
      return names_[i];
   }

   /*
   * Get nValue (number of values).
   */
   inline
   int AverageListAnalyzer::nValue(int i) const
   {  return nValue_; }

}
#endif 
