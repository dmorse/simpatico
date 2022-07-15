/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>
#include "SymmTensorAverageAnalyzer.h"
#include <ddMd/simulation/Simulation.h>
#include <util/accumulators/SymmTensorAverage.h>   
#include <util/space/Tensor.h>   
#include <util/format/Int.h>
#include <util/format/Dbl.h>
#include <util/misc/ioUtil.h>

#include <sstream>

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   SymmTensorAverageAnalyzer::SymmTensorAverageAnalyzer(Simulation& simulation) 
    : Analyzer(simulation),
      outputFile_(),
      accumulatorPtr_(0),
      nSamplePerBlock_(0),
      isInitialized_(false)
   {  setClassName("SymmTensorAverageAnalyzer"); }

   /*
   * Destructor.
   */
   SymmTensorAverageAnalyzer::~SymmTensorAverageAnalyzer() 
   {  
      if (accumulatorPtr_) {
         delete accumulatorPtr_;
      }
   }

   /*
   * Read interval and outputFileName. 
   */
   void SymmTensorAverageAnalyzer::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);
      nSamplePerBlock_ = 0;
      readOptional<int>(in,"nSamplePerBlock", nSamplePerBlock_);
      if (simulation().domain().isMaster()) {
         accumulatorPtr_ = new SymmTensorAverage;
         accumulatorPtr_->setNSamplePerBlock(nSamplePerBlock_);
      }
      isInitialized_ = true;
   }

   /*
   * Load internal state from an archive.
   */
   void SymmTensorAverageAnalyzer::loadParameters(Serializable::IArchive &ar)
   {
      loadInterval(ar);
      loadOutputFileName(ar);
      nSamplePerBlock_ = 0;
      bool isRequired = false;
      loadParameter<int>(ar, "nSamplePerBlock", nSamplePerBlock_, isRequired);
      if (simulation().domain().isMaster()) {
         accumulatorPtr_ = new SymmTensorAverage;
         ar >> *accumulatorPtr_;
         if (nSamplePerBlock_ != accumulatorPtr_->nSamplePerBlock()) {
            UTIL_THROW("Inconsistent values of nSamplePerBlock");
         }
      } else {
         accumulatorPtr_ = 0;
      }
      isInitialized_ = true;
   }

   /*
   * Save internal state to an archive.
   */
   void SymmTensorAverageAnalyzer::save(Serializable::OArchive &ar)
   {
      assert(simulation().domain().isMaster());
      assert(accumulatorPtr_);
      
      saveInterval(ar);
      saveOutputFileName(ar);
      bool isActive = (bool)nSamplePerBlock_;
      Parameter::saveOptional(ar, nSamplePerBlock_, isActive);
      ar << *accumulatorPtr_;
   }

   /*
   * Clear accumulator (do nothing on slave processors).
   */
   void SymmTensorAverageAnalyzer::clear() 
   {   
      if (simulation().domain().isMaster()){ 
         accumulatorPtr_->clear();
      }
   }
 
   /*
   * Open outputfile
   */ 
   void SymmTensorAverageAnalyzer::setup()
   {
      if (simulation().domain().isMaster()) {
         if (nSamplePerBlock_) {
            std::string filename  = outputFileName(".dat");
            simulation().fileMaster().openOutputFile(filename, outputFile_);
         }
      }
   }

   /*
   * Compute value and add to sequence.
   */
   void SymmTensorAverageAnalyzer::sample(long iStep) 
   {
      if (!isAtInterval(iStep))  {
         UTIL_THROW("Time step index is not a multiple of interval");
      }
      compute();
      if (simulation().domain().isMaster()) {
         Tensor data = value();
         accumulatorPtr_->sample(data);
         if (nSamplePerBlock_ > 0 && accumulatorPtr_->isBlockComplete()) {
            int beginStep = iStep - (nSamplePerBlock_ - 1)*interval();
            outputFile_ << Int(beginStep, 10) << "  ";
            double ave;
            int i, j;
            for (i = 0; i < Dimension; ++i) {
              for (j = 0; j <= i; ++j) {
                 ave = (*accumulatorPtr_)(i, j).blockAverage();
                 outputFile_ << Dbl(ave) << "  ";
              }
            }
            outputFile_ << "\n";
         }
      }
   }

   /*
   * Output results to file after simulation is completed.
   */
   void SymmTensorAverageAnalyzer::output()
   {
      if (simulation().domain().isMaster()) {

         // Close data (*.dat) file, if any
         if (outputFile_.is_open()) {
            outputFile_.close();
         }

         // Write parameter (*.prm) file
         FileMaster& fileMaster = simulation().fileMaster();
         fileMaster.openOutputFile(outputFileName(".prm"), outputFile_);
         ParamComposite::writeParam(outputFile_);
         outputFile_.close();

         // Write average (*.ave) file with averages for all elements
         fileMaster.openOutputFile(outputFileName(".ave"), outputFile_);
         double ave, err;
         int i, j;
         for (i = 0; i < Dimension; ++i) {
            for (j = 0; j <= i ; ++j) {
               ave = (*accumulatorPtr_)(i, j).average();
               err = (*accumulatorPtr_)(i, j).blockingError();
               outputFile_ << "Average(" << i << ", " << j << ") = ";
               outputFile_ << Dbl(ave) << " +- " << Dbl(err, 9, 2) << "\n";
            }
         }
         outputFile_.close();

         // Write average error analysis (*.aer) file
         fileMaster.openOutputFile(outputFileName(".aer"), outputFile_);
         for (i = 0; i < Dimension; ++i) {
            for (j = 0; j <= i ; ++j) {
               outputFile_ << "Element(" << i << ", " << j << "): \n\n";
               (*accumulatorPtr_)(i, j).output(outputFile_);
               outputFile_ << 
               "----------------------------------------------------------------------------\n";
            }
         }
         outputFile_.close();

      }
   }

}
