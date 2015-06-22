#ifndef DDMD_SYMM_TENSOR_AVERAGE_ANALYZER_H
#define DDMD_SYMM_TENSOR_AVERAGE_ANALYZER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/analyzers/Analyzer.h>

namespace Util { 
   class Tensor;
   class SymmTensorAverage;
}

namespace DdMd
{

   class Simulation;

   using namespace Util;

   /**
   * Analyzer that computes average of a sequence of symmetric Tensor values.
   *
   * This class evaluates the average of a sequences of symmetric Tensor values,
   * and optionally writes block averages to a data file during the run. It is
   * intended for use as a base class for classes that compute and average
   * specific symmetric-tensor-valued physical variables.
   *
   * \ingroup DdMd_Analyzer_Module
   */
   class SymmTensorAverageAnalyzer : public Analyzer
   {
   
   public:
   
      /**
      * Constructor.
      *
      * \param simulation parent Simulation object. 
      */
      SymmTensorAverageAnalyzer(Simulation& simulation);
   
      /**
      * Destructor.
      */
      virtual ~SymmTensorAverageAnalyzer(); 
   
      /**
      * Read interval, outputFileName and (optionally) nSamplePerBlock.
      *
      * The optional variable nSamplePerBlock defaults to 0, which disables
      * computation and output of block averages. Setting nSamplePerBlock = 1
      * outputs every sampled value. 

      * \param in input parameter file
      */
      virtual void readParameters(std::istream& in);
   
      /**
      * Load internal state from an input archive.
      *
      * \param ar input/loading archive
      */
      virtual void loadParameters(Serializable::IArchive &ar);

      /**
      * Save internal state to an output archive.
      *
      * \param ar output/saving archive
      */
      virtual void save(Serializable::OArchive &ar);
  
      /**
      * Clear accumulator on master, do nothing on other processors.
      */
      virtual void clear();
  
      /**
      * Setup before loop. 
      * 
      * Opens output file, if any required.
      */
      virtual void setup();

      /**
      * Compute a sampled value and add it to the accumulator.
      *
      * \param iStep MD step index
      */
      virtual void sample(long iStep);

      /**
      * Write final results to a file.
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
      virtual Tensor value() = 0;

   private:

      /// Output file stream.
      std::ofstream  outputFile_;

      /// Pointer to Average object (only instantiated on master processor)
      SymmTensorAverage *accumulatorPtr_;

      /// Number of samples per block average output.
      int nSamplePerBlock_;
   
      /// Has readParam been called?
      bool isInitialized_;
   
   };

}
#endif 
