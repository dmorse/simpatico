#ifndef SIMP_AVERAGE_MIXIN_H
#define SIMP_AVERAGE_MIXIN_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AnalyzerMixIn.h"
#include <util/archives/Serializable.h>
#include <util/global.h>

#include <iostream>
#include <string>

namespace Util { 
   class Average;
   class FileMaster;
   class ParamComposite;
}

namespace Simp
{

   using namespace Util;

   /**
   * Analyze average and block averages of a single floating point variable.
   *
   * Analyzer classes that evaluate averages in, e.g., the McMd and DdMd 
   * namespaces may be derived via multiple inheritance from both the this
   * "mixin" class and the Analyzer base class in the relevant namespace.
   * This class implements the evaluation of averages and output to file. 
   * The use of a shared MixIn avoids duplication of code in different 
   * namespaces and guarantes standardization of output formats.
   *
   * \ingroup Simp_Analysis_Module
   */
   class AverageMixIn : public AnalyzerMixIn
   {
   
   public:
   
      /**
      * Constructor.
      *
      * \param simulation  parent Simulation object. 
      */
      AverageMixIn(FileMaster& fileMaster);
   
      /**
      * Destructor.
      */
      virtual ~AverageMixIn(); 
   
      /**
      * Get value of nSamplePerBlock.
      * 
      * If nSamplePerBlock == 0, output of block averages is disabled.
      * For nSamplePerBlock > 0, the value is the number of sampled values
      * averaged in each block. 
      */
      bool nSamplePerBlock() const;

      /**
      * Does this processor have an Average accumulator?
      */
      bool hasAccumulator() const;

      /**
      * Get Average accumulator.
      *
      * \pre hasAccumulator() == true
      */
      const Average& accumulator() const;

   protected:

      /**
      * Instantiate a new Average accumulator and set nSamplePerBlock.
      *
      * \pre hasAccumulator == false
      * \pre nSamplePerBlock >= 0
      */ 
      void initializeAccumulator();

      /**
      * Clear internal state of the accumulator.
      *
      * \pre hasAccumulator == true
      */ 
      void clearAccumulator();

      /**
      * Read nSamplePerBlock parameter from file.
      *
      * \param in parameter file
      * \param composite ParamComposite that stores file format
      */ 
      void readNSamplePerBlock(std::istream& in, 
                               ParamComposite& composite);

      /**
      * Load nSamplePerBlock parameter from an archive.
      *
      * \param ar input archive
      * \param composite ParamComposite that stores file format
      */ 
      void loadNSamplePerBlock(Serializable::IArchive &ar, 
                               ParamComposite& composite);

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
      * \param interval analyzer interval, in simulation steps
      */
      void updateAccumulator(long iStep, int interval);

      /**
      * Write results of statistical analysis to files.
      *
      * \pre hasAccumulator() == true
      * \param outputFileName base output file name for analyzer
      */
      void outputAccumulator(std::string outputFileName);

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

      /// Pointer to Average object (only instantiated on master processor)
      Average *accumulatorPtr_;

      /// Number of samples per block average output.
      int nSamplePerBlock_;
   
   };

   // Inline functions

   /*
   * Does this processor have an accumulator?
   */
   inline
   bool AverageMixIn::hasAccumulator() const
   {  return (bool)(accumulatorPtr_); }

   /*
   * Get accumulator associated with a variable.
   */
   inline
   const Average& AverageMixIn::accumulator() const
   {
      UTIL_CHECK(accumulatorPtr_);
      return *accumulatorPtr_;
   }

}
#endif 
