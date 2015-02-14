/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "KineticEnergyAnalyzer.h"
#include <util/format/Int.h>
#include <util/format/Dbl.h>
#include <util/accumulators/Average.h>                    // member template 
#include <util/mpi/MpiLoader.h>
#include <util/misc/ioUtil.h>

#include <sstream>

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   KineticEnergyAnalyzer::KineticEnergyAnalyzer(Simulation& simulation) 
    : Analyzer(simulation),
      outputFile_(),
      accumulatorPtr_(0),
      nSamplePerBlock_(1),
      isInitialized_(false)
   {  setClassName("KineticEnergyAnalyzer"); }

   /*
   * Destructor.
   */
   KineticEnergyAnalyzer::~KineticEnergyAnalyzer() 
   {  
      if (accumulatorPtr_) {
         delete accumulatorPtr_;
      }
   }

   /*
   * Read interval and outputFileName. 
   */
   void KineticEnergyAnalyzer::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);
      read<int>(in,"nSamplePerBlock", nSamplePerBlock_);

      if (simulation().domain().isMaster()) {
         accumulatorPtr_ = new Average;
         accumulatorPtr_->setNSamplePerBlock(nSamplePerBlock_);
      }

      isInitialized_ = true;
   }

   /*
   * Load internal state from an archive.
   */
   void KineticEnergyAnalyzer::loadParameters(Serializable::IArchive &ar)
   {
      loadInterval(ar);
      MpiLoader<Serializable::IArchive> loader(*this, ar);
      loader.load(nSamplePerBlock_);

      if (simulation().domain().isMaster()) {
         accumulatorPtr_ = new Average;
         accumulatorPtr_->loadParameters(ar);
      }

      if (nSamplePerBlock_ != accumulatorPtr_->nSamplePerBlock()) {
         UTIL_THROW("Inconsistent values of nSamplePerBlock");
      }

      isInitialized_ = true;
   }

   /*
   * Save internal state to an archive.
   */
   void KineticEnergyAnalyzer::save(Serializable::OArchive &ar)
   {
      saveInterval(ar);
      saveOutputFileName(ar);

      if (simulation().domain().isMaster()){
         ar << *accumulatorPtr_;
      }
   }

   /*
   * Clear accumulator (do nothing on slave processors).
   */
   void KineticEnergyAnalyzer::clear() 
   {   
      if (simulation().domain().isMaster()){ 
         accumulatorPtr_->clear();
      }
   }
 
   /*
   * Compute current value.
   */
   void KineticEnergyAnalyzer::compute() 
   {  simulation().computeKineticEnergy(); }

   /*
   * Get value current value (call only on master)
   */
   double KineticEnergyAnalyzer::value() 
   {
      if (!simulation().domain().isMaster()) {
         UTIL_THROW("Error: Not master processor");
      }
      return simulation().kineticEnergy();
   }

   /*
   * Compute value.
   */
   void KineticEnergyAnalyzer::sample(long iStep) 
   {
      if (isAtInterval(iStep))  {
         compute();
         if (simulation().domain().isMaster()) {
            double data = value();
            accumulatorPtr_->sample(data);
            if (accumulatorPtr_->isBlockComplete()) {
               double block = accumulatorPtr_->blockAverage();
               outputFile_ << Int(iStep) << Dbl(block) << "\n"
            }
         }
      }
   }

   /*
   * Output results to file after simulation is completed.
   */
   void KineticEnergyAnalyzer::output()
   {
      if (simulation().domain().isMaster()) {
         simulation().fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
         ParamComposite::writeParam(outputFile_);
         outputFile_.close();

         simulation().fileMaster().openOutputFile(outputFileName(".ave"), outputFile_);
         accumulatorPtr_->output(outputFile_);
         outputFile_.close();
      }
   }

}
