#ifndef DDMD_AUTO_CORR_ANALYZER_H
#define DDMD_AUTO_CORR_ANALYZER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/analyzers/Analyzer.h>             // base class
#include <util/accumulators/AutoCorrelation.h>   // member template

#include <ddMd/simulation/Simulation.h>          // used in implementation

namespace DdMd
{

   using namespace Util;

   /**
   * Compute an autocorrelation function for a sequence of Data values.
   *
   * This template works for Data types float or double, std::complex<real>
   * with float or double real type, Util::Vector and Util::Tensor.
   *
   * \ingroup DdMd_Analyzer_Module
   */
   template <typename Data, typename Product>
   class AutoCorrAnalyzer : public Analyzer
   {
   
   public:
   
      /**
      * Constructor.
      *
      * \param simulation parent Simulation object. 
      */
      AutoCorrAnalyzer(Simulation& simulation);
   
      /**
      * Destructor.
      */
      virtual ~AutoCorrAnalyzer()
      {} 
   
      /**
      * Read interval, outputFileName and bufferCapacity from parameter file.
      *
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
      * Clear nSample counter.
      */
      virtual void clear();

      /**
      * Setup accumulator.
      */
      virtual void setup();
  
      /**
      * Compute new Data value and update accumulator.
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
      * Compute Data value, call on all processors.
      *
      * Default implementation is empty.
      */
      virtual void computeData()
      {}

      /**
      * Get current Data value, call only on master
      */
      virtual Data data() = 0;

   private:
 
      /// Output file stream
      std::ofstream  outputFile_;
      
      /// Statistical accumulator.
      AutoCorrelation<Data, Product>*  accumulatorPtr_;

      /// Buffer capacity per stage (# values stored)
      int  bufferCapacity_;

      /// Maximum stage index for descendant AutoCorrStage objects
      int  maxStageId_;

      /// Has readParam been called?
      long  isInitialized_;
   
   };

}
#endif 
