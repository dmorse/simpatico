/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AverageMixIn.h"
#include <util/accumulators/Average.h>   
#include <util/param/ParamComposite.h>   
#include <util/param/Parameter.h>   
#include <util/format/Int.h>
#include <util/format/Dbl.h>
#include <util/misc/ioUtil.h>

#include <sstream>

namespace Simp
{

   using namespace Util;

   /*
   * Constructor.
   */
   AverageMixIn::AverageMixIn(FileMaster& fileMaster) 
    : AnalyzerMixIn(fileMaster),
      accumulatorPtr_(0),
      nSamplePerBlock_(-1)
   {}

   /*
   * Destructor.
   */
   AverageMixIn::~AverageMixIn() 
   {  
      if (accumulatorPtr_) {
         delete accumulatorPtr_;
      }
   }

   /*
   * Instantiate a new Average accumulator and set nSamplerPerBlock.
   */
   void AverageMixIn::initializeAccumulator() {
      UTIL_CHECK(!hasAccumulator());
      UTIL_CHECK(nSamplePerBlock_ >= 0);
      accumulatorPtr_ = new Average;
      accumulatorPtr_->setNSamplePerBlock(nSamplePerBlock_);
   }

   /*
   * Clear accumulator.
   */
   void AverageMixIn::clearAccumulator() 
   {
      UTIL_CHECK(hasAccumulator());
      accumulatorPtr_->clear();
   }
 
   /*
   * Read nSamplePerBlock parameter from file.
   */ 
   void AverageMixIn::readNSamplePerBlock(std::istream& in,
                                          ParamComposite& composite)
   {
      // Set to zero by default
      nSamplePerBlock_ = 0;
      composite.readOptional<int>(in, "nSamplePerBlock", nSamplePerBlock_);
   }

   /*
   * Load nSamplePerBlock parameter from an archive.
   */ 
   void AverageMixIn::loadNSamplePerBlock(Serializable::IArchive &ar,
                                          ParamComposite& composite)
   {
      // Set to zero by default
      nSamplePerBlock_ = 0;
      bool isRequired = false;
      composite.loadParameter<int>(ar, "nSamplePerBlock", nSamplePerBlock_, 
                                   isRequired);
   }

   /*
   * Instantiate an Average accumulator and load data from an archive.
   */ 
   void AverageMixIn::loadAccumulator(Serializable::IArchive &ar)
   {
      UTIL_CHECK(!hasAccumulator());
      UTIL_CHECK(nSamplePerBlock_ >= 0); 
      UTIL_CHECK(accumulatorPtr_ == 0);

      accumulatorPtr_ = new Average;
      ar >> *accumulatorPtr_;
      UTIL_CHECK(accumulator().nSamplePerBlock() == nSamplePerBlock_); 
   }

   /*
   * Save nSamplePerBlock parameter to an output archive.
   */ 
   void AverageMixIn::saveNSamplePerBlock(Serializable::OArchive &ar)
   {
      UTIL_CHECK(nSamplePerBlock_ >= 0);
      bool isActive = (bool)nSamplePerBlock_;
      Parameter::saveOptional<int>(ar, nSamplePerBlock_, isActive);
   }

   /*
   * Save Average accumulator to an output archive.
   */ 
   void AverageMixIn::saveAccumulator(Serializable::OArchive &ar)
   {  ar << *accumulatorPtr_; }

   void AverageMixIn::updateAccumulator(long iStep, int interval) 
   {
      UTIL_CHECK(hasAccumulator());
      double data = value();
      accumulatorPtr_->sample(data);
      if (nSamplePerBlock_ > 0 && accumulator().isBlockComplete()) {
         UTIL_CHECK(outputFile().is_open());
         double block = accumulator().blockAverage();
         int beginStep = iStep - (nSamplePerBlock_ - 1)*interval;
         outputFile() << Int(beginStep,10) << Dbl(block, 20, 10) << "\n";
      }
   }

   /*
   * Output accumulator results to files.
   */
   void AverageMixIn::outputAccumulator(std::string outputFileName)
   {
      UTIL_CHECK(hasAccumulator());

      // Close output file.
      if (outputFile().is_open()) {
         outputFile().close();
      }

      // Write average (*.ave) file
      std::string fileName;
      fileName = outputFileName;
      fileName += ".ave"; 
      openOutputFile(fileName);
      double ave = accumulatorPtr_->average();
      double err = accumulatorPtr_->blockingError();
      outputFile() << "Average   " << Dbl(ave) << " +- " 
                  << Dbl(err, 9, 2) << "\n";
      outputFile().close();

      // Write error analysis (*.aer) file
      fileName = outputFileName;
      fileName += ".aer"; 
      openOutputFile(fileName);
      accumulatorPtr_->output(outputFile());
      outputFile().close();
   }

}
