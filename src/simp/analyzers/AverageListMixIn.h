#ifndef SIMP_AVERAGE_LIST_MIXIN_H
#define SIMP_AVERAGE_LIST_MIXIN_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <simp/analyzers/AnalyzerMixIn.h>
#include <util/accumulators/Average.h>
#include <util/containers/DArray.h>
#include <string>

namespace Util {
   class FileMaster;
   class ParamComposite;
}

namespace Simp
{

   using namespace Util;

   /**
   * Analyze averages and block averages of several real variables.
   *
   * This class evaluates the average of several sampled real variables, and
   * optionally writes block averages to a data file during a simulation. 
   * It is intended for use as a base class for Analyzers that evaluate 
   * averages and (optionally) block averages for specific physical variables.
   *
   * \ingroup Simp_Analyzer_Module
   */
   class AverageListMixIn : public AnalyzerMixIn
   {
   
   public:
   
      /**
      * Constructor.
      *
      * \param simulation  parent Simulation object. 
      */
      AverageListMixIn(FileMaster& fileMaster);
   
      /**
      * Destructor.
      */
      virtual ~AverageListMixIn(); 

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
      * \param in input parameter file
      * \param composite associated ParamComposite object
      */ 
      void readNSamplePerBlock(std::istream& in, ParamComposite& composite);

      /**
      * Load nSamplePerBlock parameter from an archive.
      *
      * \param ar input archive
      * \param composite associated ParamComposite object
      */ 
      void loadNSamplePerBlock(Serializable::IArchive &ar, ParamComposite& composite);

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
      * \param interval sampling interval
      */
      void updateAccumulators(long iStep, int interval);

      /**
      * Write final accumulator data to files.
      *  
      * \param outputFileName base output file name (without suffix)
      */
      void outputAccumulators(std::string outputFileName);

   private:

      /// Array of Average objects (only allocated on master processor)
      DArray<Average> accumulators_;

      /// Array of current values (only allocated on master processor)
      DArray<double> values_;

      /// Array of value names (only allocated on master processor)
      DArray<std::string> names_;

      /// Number of samples per block average output.
      int nSamplePerBlock_;

      /// Number of values.
      int nValue_;
 
      /// Does this processor have accumulators ?
      bool hasAccumulators_;
   
   };

   // Inline functions

   /*
   * Get nSamplePerBlock.
   */
   inline
   int AverageListMixIn::nSamplePerBlock() const
   {  
      return nSamplePerBlock_; 
   }

   /*
   * Does this processor have accumulators?
   */
   inline
   bool AverageListMixIn::hasAccumulators() const
   {  return hasAccumulators_; }

   /*
   * Get nValue (number of variables).
   */
   inline
   int AverageListMixIn::nValue() const
   {  
      UTIL_CHECK(hasAccumulators());
      return nValue_; 
   }

   /*
   * Get current value of a variable, set by compute function.
   */
   inline
   double AverageListMixIn::value(int i) const
   {
      UTIL_CHECK(hasAccumulators());
      UTIL_CHECK(i >= 0 && i < nValue_);
      return values_[i];
   }

   /*
   * Get name of specific variable.
   */
   inline
   const std::string& AverageListMixIn::name(int i) const
   {
      UTIL_CHECK(hasAccumulators());
      UTIL_CHECK(i >= 0 && i < nValue_);
      return names_[i];
   }

   /*
   * Get accumulator associated with a variable.
   */
   inline
   const Average& AverageListMixIn::accumulator(int i) const
   {
      UTIL_CHECK(hasAccumulators());
      UTIL_CHECK(i >= 0 && i < nValue_);
      return accumulators_[i];
   }

   /*
   * Set current value of a variable.
   */
   inline
   void AverageListMixIn::setValue(int i, double value)
   {
      UTIL_CHECK(hasAccumulators());
      UTIL_CHECK(i >= 0 && i < nValue_);
      values_[i] = value;
   }

}
#endif 
