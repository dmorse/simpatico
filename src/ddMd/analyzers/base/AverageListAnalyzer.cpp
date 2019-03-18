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

   // Functions declared in DdMd::Analyzer

   /*
   * Constructor.
   */
   AverageListAnalyzer::AverageListAnalyzer(Simulation& simulation) 
    : Analyzer(simulation),
      outputFile_(),
      fileMasterPtr_(&simulation.fileMaster()),
      nSamplePerBlock_(0),
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
      readNSamplePerBlock(in, *this);
      isInitialized_ = true;
   }

   /*
   * Load internal state from an archive.
   */
   void AverageListAnalyzer::loadParameters(Serializable::IArchive &ar)
   {
      loadInterval(ar);
      loadOutputFileName(ar);
      loadNSamplePerBlock(ar, *this);
      if (simulation().domain().isMaster()) {
         loadAccumulators(ar);
      }
      isInitialized_ = true;
   }

   /*
   * Save internal state to an archive.
   */
   void AverageListAnalyzer::save(Serializable::OArchive &ar)
   {
      UTIL_CHECK(simulation().domain().isMaster());
      UTIL_CHECK(hasAccumulators());
      saveInterval(ar);
      saveOutputFileName(ar);
      saveNSamplePerBlock(ar);
      saveAccumulators(ar);
   }

   /*
   * Clear accumulators (do nothing on slave processors).
   */
   void AverageListAnalyzer::clear() 
   {
      if (hasAccumulators()) {
         clearAccumulators();
      }
   }

   /*
   * Setup before simulation.
   */ 
   void AverageListAnalyzer::setup()
   {
      if (hasAccumulators()) {
         if (nSamplePerBlock()) {
            openOutputFile(outputFileName(".dat"));
         }
      }
   }

   /*
   * Compute and sample current values.
   */
   void AverageListAnalyzer::sample(long iStep) 
   {
      if (!isAtInterval(iStep)) {
         UTIL_THROW("Time step index is not a multiple of interval");
      }
      compute();
      if (hasAccumulators()) {
         updateAccumulators(iStep);
      }
   }

   /*
   * Output results after a simulation is completed.
   */
   void AverageListAnalyzer::output()
   {
      if (hasAccumulators()) {

         // Close data (*.dat) file, if any
         if (outputFile_.is_open()) {
            outputFile_.close();
         }
   
         // Write parameter (*.prm) file
         openOutputFile(outputFileName(".prm"));
         ParamComposite::writeParam(outputFile_);
         outputFile_.close();

         // Write average (*.ave) and error analysis (*.aer) files
         outputAccumulators();

      }
   }

   // Utility functions (migrate to base class)

   /**
   * Set nValue and allocate arrays with dimensions nValue.
   */ 
   void AverageListAnalyzer::initializeAccumulators(int nValue) 
   {
      UTIL_CHECK(nValue > 0);
      UTIL_CHECK(nValue_ == 0);
      accumulators_.allocate(nValue);
      names_.allocate(nValue);
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
   * Clear accumulators.
   */
   void AverageListAnalyzer::clearAccumulators() 
   {
      UTIL_CHECK(hasAccumulators());
      UTIL_CHECK(nValue_ > 0);
      for (int i = 0; i < nValue_; ++i) {
         accumulators_[i].clear();
      }
   }

   /*
   * Read nSamplePerBlock parameter from file.
   */ 
   void AverageListAnalyzer::readNSamplePerBlock(std::istream& in,
                                                 ParamComposite& composite)
   {
      nSamplePerBlock_ = 0;
      composite.readOptional<int>(in, "nSamplePerBlock", nSamplePerBlock_);
   }

   /*
   * Load nSamplePerBlock parameter from an archive.
   */ 
   void AverageListAnalyzer::loadNSamplePerBlock(Serializable::IArchive &ar,
                                                 ParamComposite& composite)
   {
      nSamplePerBlock_ = 0;
      bool isRequired = false;
      composite.loadParameter<int>(ar, "nSamplePerBlock", nSamplePerBlock_, 
                                   isRequired);
   }

   /*
   * Instantiate an Average accumulator and load data from an archive.
   */ 
   void AverageListAnalyzer::loadAccumulators(Serializable::IArchive &ar)
   {
      UTIL_CHECK(nSamplePerBlock_ > 0);
      UTIL_CHECK(!hasAccumulators());
      int nValue;
      ar >> nValue;
      UTIL_CHECK(nValue > 0);
      initializeAccumulators(nValue);
      UTIL_CHECK(nValue_ == nValue);
      ar >> accumulators_;
      ar >> names_;
      for (int i = 0; i < nValue_; ++i) {
         UTIL_CHECK(accumulator(i).nSamplePerBlock() == nSamplePerBlock_);
      }
   }

   /*
   * Save nSamplePerBlock parameter to an output archive.
   */ 
   void AverageListAnalyzer::saveNSamplePerBlock(Serializable::OArchive &ar)
   {
      bool isActive = (bool)nSamplePerBlock_;
      Parameter::saveOptional<int>(ar, nSamplePerBlock_, isActive);
   }

   /*
   * Save Average accumulator to an output archive.
   */ 
   void AverageListAnalyzer::saveAccumulators(Serializable::OArchive &ar)
   {  
      UTIL_CHECK(hasAccumulators()); 
      UTIL_CHECK(nValue_ > 0);
      ar << nValue_;
      ar << accumulators_;
      ar << names_;
   }

   /*
   * Update accumulators for all current values.
   */
   void AverageListAnalyzer::updateAccumulators(long iStep) 
   {
      UTIL_CHECK(hasAccumulators());

      // Decide whether to output block averages
      bool outputBlocks = false;
      if (nSamplePerBlock_ > 0) { 
         if (accumulators_[0].isBlockComplete() > 0) {
            outputBlocks = true;
         }
      }

      // If outputBlocks, write step counter
      if (outputBlocks) {
         UTIL_CHECK(outputFile_.is_open());
         int beginStep = iStep - (nSamplePerBlock_ - 1)*interval();
         outputFile_ << Int(beginStep);
      }

      // Loop over values
      for (int i = 0; i < nValue(); ++i) {
         double data = value(i);
         accumulators_[i].sample(data);
         if (outputBlocks) {
            UTIL_CHECK(accumulators_[i].isBlockComplete());
            double block = accumulators_[i].blockAverage();
            outputFile_ << Dbl(block);
         }
      }

      if (outputBlocks) {
         outputFile_ << "\n";
      }
   }

   /*
   * Output results to file after simulation is completed.
   */
   void AverageListAnalyzer::outputAccumulators()
   {
      UTIL_CHECK(hasAccumulators());

      // Close data (*.dat) file, if any
      if (outputFile_.is_open()) {
         outputFile_.close();
      }

      // Write average (*.ave) file
      openOutputFile(outputFileName(".ave"));
      double ave, err;
      for (int i = 0; i < nValue_; ++i) {
         ave = accumulators_[i].average();
         err = accumulators_[i].blockingError();
         outputFile_ << names_[i] << "   ";
         outputFile_ << Dbl(ave) << " +- " << Dbl(err, 9, 2) << "\n";
      }
      outputFile_.close();

      // Write error analysis (*.aer) file
      openOutputFile(outputFileName(".aer"));
      for (int i = 0; i < nValue_; ++i) {
         accumulators_[i].output(outputFile_);
      }
      outputFile_.close();
   }

   /*
   * Open the ouput file with specified base name.
   */
   void 
   AverageListAnalyzer::openOutputFile(std::string name)
   {  fileMasterPtr_->openOutputFile(name, outputFile_); }

}
