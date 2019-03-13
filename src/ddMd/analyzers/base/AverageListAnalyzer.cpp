/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AverageListAnalyzer.h"
#include <ddMd/simulation/Simulation.h>
#include <util/accumulators/Average.h>   
#include <util/format/Int.h>
#include <util/format/Dbl.h>
#include <util/misc/FileMaster.h>
#include <util/misc/ioUtil.h>

#include <sstream>

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   AverageListAnalyzer::AverageListAnalyzer(Simulation& simulation) 
    : Analyzer(simulation),
      nSamplePerBlock_(0),
      outputFile_(),
      fileMasterPtr_(&simulation.fileMaster()),
      nValue_(0),
      hasAccumulators_(false),
      isInitialized_(false)
   {  setClassName("AverageListAnalyzer"); }

   /*
   * Destructor.
   */
   AverageListAnalyzer::~AverageListAnalyzer() 
   {}

   /*
   * Read interval and outputFileName. 
   */
   void AverageListAnalyzer::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);
      nSamplePerBlock_ = 0;
      readOptional<int>(in, "nSamplePerBlock", nSamplePerBlock_);
      isInitialized_ = true;
   }

   /*
   * Load internal state from an archive.
   */
   void AverageListAnalyzer::loadParameters(Serializable::IArchive &ar)
   {
      loadInterval(ar);
      loadOutputFileName(ar);
      nSamplePerBlock_ = 0;
      bool isRequired = false;
      loadParameter<int>(ar, "nSamplePerBlock", nSamplePerBlock_, 
                         isRequired);
      if (hasAccumulators_) {
         ar >> nValue_;
         if (nValue_ > 0) {
            ar >> accumulators_;
            ar >> names_;
         }
         for (int i = 0; i < nValue_; ++i) {
            UTIL_CHECK(accumulators_[i].nSamplePerBlock() == nSamplePerBlock_);
         }
      }
      isInitialized_ = true;
   }

   /*
   * Save internal state to an archive.
   */
   void AverageListAnalyzer::save(Serializable::OArchive &ar)
   {
      UTIL_CHECK(hasAccumulators_);

      saveInterval(ar);
      saveOutputFileName(ar);
      bool isActive = (bool)nSamplePerBlock_;
      Parameter::saveOptional<int>(ar, nSamplePerBlock_, isActive);
      if (hasAccumulators_ > 0) {
         UTIL_CHECK(nValue_ > 0);
         ar << nValue_;
         ar << accumulators_;
         ar << names_;
      }
   }

   /*
   * Clear accumulator (do nothing on slave processors).
   */
   void AverageListAnalyzer::clear() 
   {
      if (hasAccumulators_ > 0) {
         for (int i = 0; i < nValue_; ++i) {
            accumulators_[i].clear();
         }
      }
   }

   /**
   * Set nValue and allocate arrays with dimensions nValue.
   */ 
   void AverageListAnalyzer::setNValue(int nValue) 
   {
      UTIL_CHECK(nValue > 0);
      UTIL_CHECK(nValue_ == 0);
      names_.allocate(nValue);
      accumulators_.allocate(nValue);
      nValue_ = nValue;
      hasAccumulators_ = true;
   }

   void AverageListAnalyzer::setName(int i, std::string name) 
   {
      UTIL_CHECK(hasAccumulators_ > 0);
      UTIL_CHECK(i > 0 && i < nValue_);
      names_[i] = name;
   }

   /*
   * Open outputfile
   */ 
   void AverageListAnalyzer::setup()
   {
      UTIL_CHECK(nSamplePerBlock_ >= 0);
      if (hasAccumulators_) {
         if (nSamplePerBlock_) {
            openOutputFile(outputFileName(".dat"), outputFile_);
         }
      }
   }

   /*
   * Compute value.
   */
   void AverageListAnalyzer::sample(long iStep) 
   {
      if (!isAtInterval(iStep)) {
         UTIL_THROW("Time step index is not a multiple of interval");
      }
      compute();
      if (hasAccumulators_) {
         if (nSamplePerBlock_ > 0 
                                && accumulators_[0].isBlockComplete() > 0) {
            int beginStep = iStep - (nSamplePerBlock_ - 1)*interval();
            outputFile_ << Int(beginStep);
         }
         for (int i = 0; i < nValue(); ++i) {
            double data = value(i);
            accumulators_[i].sample(data);
            if (nSamplePerBlock_ > 0 && accumulators_[i].isBlockComplete()){
               double block = accumulators_[i].blockAverage();
               outputFile_ << Dbl(block);
            }
         }
         outputFile_ << "\n";
      }
   }

   /*
   * Output results to file after simulation is completed.
   */
   void AverageListAnalyzer::output()
   {
      if (hasAccumulators_) {
         // Close data (*.dat) file, if any
         if (outputFile_.is_open()) {
            outputFile_.close();
         }

         // Write parameter (*.prm) file
         openOutputFile(outputFileName(".prm"), outputFile_);
         ParamComposite::writeParam(outputFile_);
         outputFile_.close();

         // Write average (*.ave) file
         openOutputFile(outputFileName(".ave"), outputFile_);
         double ave, err;
         for (int i = 0; i < nValue_; ++i) {
            ave = accumulators_[i].average();
            err = accumulators_[i].blockingError();
            outputFile_ << names_[i] << "   ";
            outputFile_ << Dbl(ave) << " +- " << Dbl(err, 9, 2) << "\n";
         }
         outputFile_.close();

         // Write error analysis (*.aer) file
         openOutputFile(outputFileName(".aer"), outputFile_);
         for (int i = 0; i < nValue_; ++i) {
            accumulators_[i].output(outputFile_);
         }
         outputFile_.close();
      }
   }

   void 
   AverageListAnalyzer::openOutputFile(std::string name, 
                                       std::ofstream& file) 
   {  fileMasterPtr_->openOutputFile(name, outputFile_); }

}
