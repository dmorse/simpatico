/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AverageListMixIn.h"
#include <util/format/Int.h>
#include <util/format/Dbl.h>
#include <util/misc/FileMaster.h>
#include <util/misc/ioUtil.h>

#include <sstream>

namespace Simp
{

   using namespace Util;

   /*
   * Constructor.
   */
   AverageListMixIn::AverageListMixIn(FileMaster& fileMaster) 
    : AnalyzerMixIn(fileMaster),
      nSamplePerBlock_(0),
      nValue_(0),
      hasAccumulators_(false)
   {}

   /*
   * Destructor.
   */
   AverageListMixIn::~AverageListMixIn() 
   {}

   /**
   * Set nValue and allocate arrays with dimensions nValue.
   */ 
   void AverageListMixIn::initializeAccumulators(int nValue) 
   {
      UTIL_CHECK(nValue > 0);
      UTIL_CHECK(nValue_ == 0);
      UTIL_CHECK(nSamplePerBlock_ >= 0);
      accumulators_.allocate(nValue);
      names_.allocate(nValue);
      values_.allocate(nValue);
      nValue_ = nValue;
      hasAccumulators_ = true;
      for (int i = 0; i < nValue_; ++i) {
         accumulators_[i].setNSamplePerBlock(nSamplePerBlock_);
      }
      clearAccumulators();
   }

   void AverageListMixIn::setName(int i, std::string name) 
   {
      UTIL_CHECK(hasAccumulators_);
      UTIL_CHECK(i >= 0 && i < nValue_);
      names_[i] = name;
   }

   /*
   * Clear accumulators.
   */
   void AverageListMixIn::clearAccumulators() 
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
   void AverageListMixIn::readNSamplePerBlock(std::istream& in,
                                              ParamComposite& composite)
   {
      nSamplePerBlock_ = 0;
      composite.readOptional<int>(in, "nSamplePerBlock", nSamplePerBlock_);
   }

   /*
   * Load nSamplePerBlock parameter from an archive.
   */ 
   void AverageListMixIn::loadNSamplePerBlock(Serializable::IArchive &ar,
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
   void AverageListMixIn::loadAccumulators(Serializable::IArchive &ar)
   {
      UTIL_CHECK(nSamplePerBlock_ >= 0);
      UTIL_CHECK(!hasAccumulators());
      int nValue;
      ar >> nValue;
      UTIL_CHECK(nValue > 0);
      initializeAccumulators(nValue);
      UTIL_CHECK(nValue_ == nValue);
      UTIL_CHECK(accumulators_.capacity() == nValue_);
      UTIL_CHECK(names_.capacity() == nValue_);
      UTIL_CHECK(values_.capacity() == nValue_);
      ar >> accumulators_;
      ar >> names_;
      for (int i = 0; i < nValue_; ++i) {
         UTIL_CHECK(accumulator(i).nSamplePerBlock() == nSamplePerBlock_);
      }
   }

   /*
   * Save nSamplePerBlock parameter to an output archive.
   */ 
   void AverageListMixIn::saveNSamplePerBlock(Serializable::OArchive &ar)
   {
      bool isActive = (bool)nSamplePerBlock_;
      Parameter::saveOptional<int>(ar, nSamplePerBlock_, isActive);
   }

   /*
   * Save Average accumulator to an output archive.
   */ 
   void AverageListMixIn::saveAccumulators(Serializable::OArchive &ar)
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
   void AverageListMixIn::updateAccumulators(long iStep, int interval) 
   {
      UTIL_CHECK(hasAccumulators());
      UTIL_CHECK(accumulators_.capacity() == nValue_);

      // Update accumulators.
      for (int i = 0; i < nValue(); ++i) {
         double data = value(i);
         accumulators_[i].sample(data);
      }

      // Output block averages
      if (nSamplePerBlock_ > 0) { 
         if (accumulators_[0].isBlockComplete()) {
            UTIL_CHECK(outputFile().is_open());
            int beginStep = iStep - (nSamplePerBlock_ - 1)*interval;
            outputFile() << Int(beginStep);
            for (int i = 0; i < nValue(); ++i) {
               UTIL_CHECK(accumulators_[i].isBlockComplete());
               double block = accumulators_[i].blockAverage();
               outputFile() << Dbl(block);
            }
            outputFile() << "\n";
         }
      }

   }

   /*
   * Output results to file after simulation is completed.
   */
   void AverageListMixIn::outputAccumulators(std::string outputFileName)
   {
      UTIL_CHECK(hasAccumulators());

      // Close data (*.dat) file, if any
      if (outputFile().is_open()) {
         outputFile().close();
      }

      // Compute maximum length of name field
      int nameWidth = 0;
      int length;
      for (int i = 0; i < nValue_; ++i) {
         length = names_[i].length();
         if (length > nameWidth) {
            nameWidth = length;
         }
      }
      nameWidth += 2;

      // Write average (*.ave) file
      std::string fileName;
      fileName = outputFileName;
      fileName += ".ave";
      openOutputFile(fileName);
      double ave, err;
      for (int i = 0; i < nValue_; ++i) {
         ave = accumulators_[i].average();
         err = accumulators_[i].blockingError();
         outputFile() << " " << std::left << std::setw(nameWidth) 
                      << names_[i] << "   ";
         outputFile() << Dbl(ave) << " +- " << Dbl(err, 9, 2) << "\n";
      }
      outputFile().close();

      // Write error analysis (*.aer) file
      fileName = outputFileName;
      fileName += ".aer";
      openOutputFile(fileName);
      std::string line; 
      line = 
      "---------------------------------------------------------------------";
      for (int i = 0; i < nValue_; ++i) {
         outputFile() << line << std::endl;
         outputFile() << names_[i] << " :" << std::endl;
         accumulators_[i].output(outputFile());
         outputFile() << std::endl;
      }
      outputFile().close();

      // Write data format file (*.dfm) file
      fileName = outputFileName;
      fileName += ".dfm";
      openOutputFile(fileName);
      outputFile() << "Value = " << nValue() << std::endl;
      outputFile() << "iStep  ";
      for (int i = 0; i < nValue_; ++i) {
         outputFile() << names_[i] << "  ";
      }
      outputFile() << std::endl;
      outputFile().close();

   }

}
