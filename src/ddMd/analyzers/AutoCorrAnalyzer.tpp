/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AutoCorrAnalyzer.h"
#include <ddMd/simulation/Simulation.h>
#include <util/accumulators/AutoCorrelation.tpp>
#include <util/misc/ioUtil.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

#include <sstream>

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   template <typename Data, typename Product>
   AutoCorrAnalyzer<Data, Product>::AutoCorrAnalyzer(Simulation& simulation) 
    : Analyzer(simulation),
      accumulatorPtr_(0),
      bufferCapacity_(-1),
      maxStageId_(10),
      isInitialized_(false)
   {  setClassName("AutoCorrAnalyzer"); }

   /*
   * Read interval, outputFileName and bufferCapacity_.
   */
   template <typename Data, typename Product>
   void AutoCorrAnalyzer<Data, Product>::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);
      read<int>(in,"bufferCapacity", bufferCapacity_);
      if (simulation().domain().isMaster()) {
         accumulatorPtr_ = new AutoCorrelation<Tensor, double>;
         accumulatorPtr_->setParam(bufferCapacity_, maxStageId_);
      }
      isInitialized_ = true;
   }


   /*
   * Load internal state from an archive.
   */
   template <typename Data, typename Product>
   void AutoCorrAnalyzer<Data, Product>::loadParameters(Serializable::IArchive &ar)
   {
      loadInterval(ar);
      loadOutputFileName(ar);
      loadParameter(ar, "bufferCapacity", bufferCapacity_);

      if (simulation().domain().isMaster()) {
         accumulatorPtr_ = new AutoCorrelation<Tensor, double>;
         ar >> *accumulatorPtr_;
      }

      isInitialized_ = true;
   }

   /*
   * Save internal state to an archive.
   */
   template <typename Data, typename Product>
   void AutoCorrAnalyzer<Data, Product>::save(Serializable::OArchive &ar)
   {
      saveInterval(ar);
      saveOutputFileName(ar);
      ar & bufferCapacity_;
      if (simulation().domain().isMaster()) {
         if (!accumulatorPtr_) {
            UTIL_THROW("Null accumulatorPtr_ on master");
         }
         ar << *accumulatorPtr_;
      }
   }
  
   /*
   * Clear accumulator.
   */
   template <typename Data, typename Product>
   void AutoCorrAnalyzer<Data, Product>::clear() 
   {
      if (!isInitialized_) {
         UTIL_THROW("Error: Object not initialized");
      }
      if (simulation().domain().isMaster()) {
         if (!accumulatorPtr_) {
            UTIL_THROW("Null accumulatorPtr_ on master");
         }
         accumulatorPtr_->clear();  
      }
   }

   /*
   * Setup before simulation.
   */
   template <typename Data, typename Product>
   void AutoCorrAnalyzer<Data, Product>::setup()
   {
      if (!isInitialized_) {
         UTIL_THROW("Object not initialized.");
      }
      // clear();
   }

   /*
   * Sample one Data value.
   */
   template <typename Data, typename Product>
   void AutoCorrAnalyzer<Data, Product>::sample(long iStep) 
   {  
      if (!isAtInterval(iStep))  {
         UTIL_THROW("Time step index is not a multiple of interval");
      }
      computeData();
      if (simulation().domain().isMaster()) {
         if (!accumulatorPtr_) {
            UTIL_THROW("Null accumulatorPtr_ on master");
         }
         Data value = data();
         accumulatorPtr_->sample(value);
      }
   }

   /*
   * Output autocorrelation function to file.
   */
   template <typename Data, typename Product>
   void AutoCorrAnalyzer<Data, Product>::output() 
   {
      if (simulation().domain().isMaster()) {
         if (!accumulatorPtr_) {
            UTIL_THROW("Null accumulatorPtr_ on master");
         }

         simulation().fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
         writeParam(outputFile_);
         outputFile_ << std::endl;
         outputFile_ << "bufferCapacity  " << accumulatorPtr_->bufferCapacity() << std::endl;
         outputFile_ << "nSample         " << accumulatorPtr_->nSample() << std::endl;
         outputFile_ << std::endl;
         outputFile_ << "Format of *.dat file" << std::endl;
         outputFile_ << "[int time delay (samples)]  [double autocorrelation function]"
                     << std::endl;
         outputFile_ << std::endl;
         outputFile_.close();

         // Write xy autocorrelation function to data file
         simulation().fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
         accumulatorPtr_->output(outputFile_);
         outputFile_.close();
      }
   }

}
