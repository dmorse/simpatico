#ifndef DDMD_TENSOR_AVERAGE_ANALYZER_H
#define DDMD_TENSOR_AVERAGE_ANALYZER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/analyzers/Analyzer.h>  // Base class header

namespace Util { 
   class Tensor;
   class TensorAverage;
}

namespace DdMd
{

   class Simulation;

   using namespace Util;

   /**
   * Analyzer that computes averages and block averages of a single Tensor.
   *
   * This class evaluates the average of a sequence of Tensor values, and
   * optionally writes block averages to a data file during the run. It is
   * intended for use as a base class for classes that evaluate averages 
   * for specific tensor-valued physical variables.
   *
   * \ingroup DdMd_Analyzer_Module
   */
   class TensorAverageAnalyzer : public Analyzer
   {
   
   public:
   
      /**
      * Constructor.
      *
      * \param simulation parent Simulation object. 
      */
      TensorAverageAnalyzer(Simulation& simulation);
   
      /**
      * Destructor.
      */
      virtual ~TensorAverageAnalyzer(); 
   
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
      * Load internal state from an archive.
      *
      * \param ar  input/loading archive
      */
      virtual void loadParameters(Serializable::IArchive &ar);

      /**
      * Save internal state to an archive.
      *
      * \param ar  output/saving archive
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
      * Compute a sampled value and add update accumulator.
      *
      * If block averaging is enabled (nSamplePerBlock > 0), then
      * this function also periodically outputs block averages.
      *
      * \param iStep  MD time step index
      */
      virtual void sample(long iStep);

      /**
      * Write final average and error analysis to file.
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
      TensorAverage *accumulatorPtr_;

      /// Number of samples per block average output.
      int nSamplePerBlock_;
   
      /// Has readParam been called?
      bool isInitialized_;
   
   };

}
#endif 
