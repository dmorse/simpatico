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
#include <util/misc/FileMaster.h>
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
      fileMasterPtr_(&simulation.fileMaster()),
      accumulatorPtr_(0),
      nSamplePerBlock_(0)
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
      readNSamplePerBlock(in, *this);
      if (simulation().domain().isMaster()) {
         initializeAccumulator();
      }
   }

   /*
   * Load internal state from an archive.
   */
   void AverageAnalyzer::loadParameters(Serializable::IArchive &ar)
   {
      loadInterval(ar);
      loadOutputFileName(ar);
      loadNSamplePerBlock(ar, *this);
      if (simulation().domain().isMaster()) {
         loadAccumulator(ar);
      }
   }

   /*
   * Save internal state to an archive.
   */
   void AverageAnalyzer::save(Serializable::OArchive &ar)
   {
      UTIL_CHECK(simulation().domain().isMaster());
      UTIL_CHECK(hasAccumulator());

      saveInterval(ar);
      saveOutputFileName(ar);
      saveNSamplePerBlock(ar);
      saveAccumulator(ar);
   }

   /*
   * Clear accumulator (do nothing on slave processors).
   */
   void AverageAnalyzer::clear() 
   {   
      if (simulation().domain().isMaster()){ 
         clearAccumulator();
      }
   }
 
   /*
   * Open outputfile
   */ 
   void AverageAnalyzer::setup()
   {
      if (simulation().domain().isMaster()) {
         UTIL_CHECK(hasAccumulator());
         if (nSamplePerBlock() > 0) {
            openOutputFile(outputFileName(".dat"));
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
         updateAccumulator(iStep);
      }
   }

   /*
   * Output results to several files after simulation is completed.
   */
   void AverageAnalyzer::output()
   {
      if (hasAccumulator()) {

         // Close data (*.dat) file, if any
         if (outputFile().is_open()) {
            outputFile().close();
         }

         // Write parameter (*.prm) file
         openOutputFile(outputFileName(".prm"));
         ParamComposite::writeParam(outputFile());
         outputFile().close();

         // Write average (*.ave) and error analysis (*.aer) files
         outputAccumulator();

      }
   }

   // Utility functions (Migrate to base class)

   /*
   * Instantiate a new Average accumulator and set nSamplerPerBlock.
   */
   void AverageAnalyzer::initializeAccumulator() {
      UTIL_CHECK(!hasAccumulator());
      UTIL_CHECK(nSamplePerBlock_ > 0);
      accumulatorPtr_ = new Average;
      accumulatorPtr_->setNSamplePerBlock(nSamplePerBlock_);
   }

   /*
   * Clear accumulator.
   */
   void AverageAnalyzer::clearAccumulator() 
   {
      UTIL_CHECK(hasAccumulator());
      accumulatorPtr_->clear();
   }
 
   /*
   * Read nSamplePerBlock parameter from file.
   */ 
   void AverageAnalyzer::readNSamplePerBlock(std::istream& in, 
                                             ParamComposite& composite)
   {
      nSamplePerBlock_ = 0;
      composite.readOptional<int>(in, "nSamplePerBlock", nSamplePerBlock_);
   }

   /*
   * Load nSamplePerBlock parameter from an archive.
   */ 
   void AverageAnalyzer::loadNSamplePerBlock(Serializable::IArchive &ar, 
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
   void AverageAnalyzer::loadAccumulator(Serializable::IArchive &ar)
   {
      UTIL_CHECK(!hasAccumulator());
      UTIL_CHECK(nSamplePerBlock_ > 0); 
      accumulatorPtr_ = new Average;
      ar >> *accumulatorPtr_;
      UTIL_CHECK(nSamplePerBlock_ != accumulator().nSamplePerBlock()); 
   }

   /*
   * Save nSamplePerBlock parameter to an output archive.
   */ 
   void AverageAnalyzer::saveNSamplePerBlock(Serializable::OArchive &ar)
   {
      bool isActive = (bool)nSamplePerBlock_;
      Parameter::saveOptional<int>(ar, nSamplePerBlock_, isActive);
   }

   /*
   * Save Average accumulator to an output archive.
   */ 
   void AverageAnalyzer::saveAccumulator(Serializable::OArchive &ar)
   {  ar << *accumulatorPtr_; }

   void AverageAnalyzer::updateAccumulator(long iStep) 
   {
      UTIL_CHECK(hasAccumulator());
      double data = value();
      accumulatorPtr_->sample(data);
      if (nSamplePerBlock_ > 0 && accumulator().isBlockComplete()) {
         UTIL_CHECK(outputFile().is_open());
         double block = accumulator().blockAverage();
         int beginStep = iStep - (nSamplePerBlock_ - 1)*interval();
         outputFile() << Int(beginStep,10) << Dbl(block, 20, 10) << "\n";
      }
   }

   /*
   * Output accumulator results to files.
   */
   void AverageAnalyzer::outputAccumulator()
   {
      UTIL_CHECK(hasAccumulator());

      // Close output file.
      if (outputFile().is_open()) {
         outputFile().close();
      }

      // Write average (*.ave) file
      openOutputFile(outputFileName(".ave"));
      double ave = accumulatorPtr_->average();
      double err = accumulatorPtr_->blockingError();
      outputFile() << "Average   " << Dbl(ave) << " +- " 
                  << Dbl(err, 9, 2) << "\n";
      outputFile().close();

      // Write error analysis (*.aer) file
      openOutputFile(outputFileName(".aer"));
      accumulatorPtr_->output(outputFile());
      outputFile().close();
   }

   /*
   * Open the output file.
   */
   void AverageAnalyzer::openOutputFile(std::string name)
   {  fileMasterPtr_->openOutputFile(name, outputFile_); }

}
