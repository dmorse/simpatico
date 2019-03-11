/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AverageAnalyzer.h"
#include <ddMd/simulation/Simulation.h>
#include <util/accumulators/Average.h>   
#include <util/format/Int.h>
#include <util/format/Dbl.h>
#include <util/mpi/MpiLoader.h>
#include <util/misc/ioUtil.h>

#include <sstream>

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   AverageAnalyzer::AverageAnalyzer(Simulation& simulation) 
    : Analyzer(simulation),
      outputFile_(),
      accumulatorPtr_(0),
      nSamplePerBlock_(0),
      isInitialized_(false)
   {  setClassName("AverageAnalyzer"); }

   /*
   * Destructor.
   */
   AverageAnalyzer::~AverageAnalyzer() 
   {  
      if (accumulatorPtr_) {
         delete accumulatorPtr_;
      }
   }

   /*
   * Read interval and outputFileName. 
   */
   void AverageAnalyzer::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);
      nSamplePerBlock_ = 0;
      readOptional<int>(in, "nSamplePerBlock", nSamplePerBlock_);

      if (simulation().domain().isMaster()) {
         accumulatorPtr_ = new Average;
         accumulatorPtr_->setNSamplePerBlock(nSamplePerBlock_);
      }

      isInitialized_ = true;
   }

   /*
   * Load internal state from an archive.
   */
   void AverageAnalyzer::loadParameters(Serializable::IArchive &ar)
   {
      loadInterval(ar);
      loadOutputFileName(ar);
      nSamplePerBlock_ = 0;
      bool isRequired = false;
      loadParameter<int>(ar, "nSamplePerBlock", nSamplePerBlock_, 
                         isRequired);

      if (simulation().domain().isMaster()) {
         accumulatorPtr_ = new Average;
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
   void AverageAnalyzer::save(Serializable::OArchive &ar)
   {
      assert(simulation().domain().isMaster());
      assert(accumulatorPtr_);

      saveInterval(ar);
      saveOutputFileName(ar);
      bool isActive = (bool)nSamplePerBlock_;
      Parameter::saveOptional<int>(ar, nSamplePerBlock_, isActive);
      ar << *accumulatorPtr_;
   }

   /*
   * Clear accumulator (do nothing on slave processors).
   */
   void AverageAnalyzer::clear() 
   {   
      if (simulation().domain().isMaster()){ 
         accumulatorPtr_->clear();
      }
   }
 
   /*
   * Open outputfile
   */ 
   void AverageAnalyzer::setup()
   {
      if (simulation().domain().isMaster()) {
         if (nSamplePerBlock_) {
            std::string filename  = outputFileName(".dat");
            simulation().fileMaster().openOutputFile(filename, outputFile_);
         }
      }
   }

   /*
   * Compute value.
   */
   void AverageAnalyzer::sample(long iStep) 
   {
      if (!isAtInterval(iStep)) {
         UTIL_THROW("Time step index is not a multiple of interval");
      }
      compute();
      if (simulation().domain().isMaster()) {
         double data = value();
         accumulatorPtr_->sample(data);
         if (nSamplePerBlock_ > 0 && accumulatorPtr_->isBlockComplete()) {
            double block = accumulatorPtr_->blockAverage();
            int beginStep = iStep - (nSamplePerBlock_ - 1)*interval();
            outputFile_ << Int(beginStep) << Dbl(block) << "\n";
         }
      }
   }

   /*
   * Output results to file after simulation is completed.
   */
   void AverageAnalyzer::output()
   {
      if (simulation().domain().isMaster()) {
         // Close data (*.dat) file, if any
         if (outputFile_.is_open()) {
            outputFile_.close();
         }
         // Write parameter (*.prm) file
         simulation().fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
         ParamComposite::writeParam(outputFile_);
         outputFile_.close();
         // Write average (*.ave) file
         simulation().fileMaster().openOutputFile(outputFileName(".ave"), outputFile_);
         double ave = accumulatorPtr_->average();
         double err = accumulatorPtr_->blockingError();
         outputFile_ << "Average   " << Dbl(ave) << " +- " << Dbl(err, 9, 2) << "\n";
         outputFile_.close();
         // Write error analysis (*.aer) file
         simulation().fileMaster().openOutputFile(outputFileName(".aer"), outputFile_);
         accumulatorPtr_->output(outputFile_);
         outputFile_.close();
      }
   }

}
