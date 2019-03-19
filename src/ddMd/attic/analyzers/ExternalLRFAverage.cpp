/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ExternalLRFAverage.h"
#include <ddMd/potentials/pair/PairPotential.h>
#include <ddMd/potentials/external/ExternalPotential.h>
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
   ExternalLRFAverage::ExternalLRFAverage(Simulation& simulation) 
    : Analyzer(simulation),
      outputFile_(),
      accumulator_(NULL),
      nSamplePerBlock_(1),
      isInitialized_(false)
   {  setClassName("ExternalLRFAverage"); }

   /*
   * Destructor.
   */
   ExternalLRFAverage::~ExternalLRFAverage() 
   {  
      if(accumulator_ != NULL) {
         delete accumulator_;
      }
   }

   /*
   * Read interval and outputFileName. 
   */
   void ExternalLRFAverage::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);
      read<int>(in,"nSamplePerBlock", nSamplePerBlock_);

      if (simulation().domain().isMaster()) {
         accumulator_ = new Average;
         accumulator_->setNSamplePerBlock(nSamplePerBlock_);
      }

      isInitialized_ = true;
   }

   /*
   * Load internal state from an archive.
   */
   void ExternalLRFAverage::loadParameters(Serializable::IArchive &ar)
   {
      loadInterval(ar);
      loadOutputFileName(ar);
      loadParameter<int>(ar,"nSamplePerBlock", nSamplePerBlock_);
      if (simulation().domain().isMaster()) {
         accumulator_ = new Average;
         accumulator_->loadParameters(ar);
         if (nSamplePerBlock_ != accumulator_->nSamplePerBlock()) {
            UTIL_THROW("Inconsistent values of nSamplePerBlock");
         }
      }
      isInitialized_ = true;
   }

   /*
   * Save internal state to an archive.
   */
   void ExternalLRFAverage::save(Serializable::OArchive &ar)
   {
      assert(simulation().domain().isMaster());
      saveInterval(ar);
      saveOutputFileName(ar);
      ar << nSamplePerBlock_;
      ar << *accumulator_;
   }

   /*
   * Reset nSample.
   */
   void ExternalLRFAverage::clear() 
   {   
      if (simulation().domain().isMaster()){ 
         accumulator_->clear();
      }
   }

   /*
   * Dump configuration to file
   */
   void ExternalLRFAverage::sample(long iStep) 
   {
      if (isAtInterval(iStep))  {
         double lrf = 0.0;
         double ext = 0.0;
         simulation().computePotentialEnergies();
         #ifdef SIMP_EXTERNAL
         ext = simulation().externalPotential().get("externalParameter");
         simulation().externalPotential().set("externalParameter",1.0);
         if (simulation().hasExternal()) {
            lrf =  simulation().externalPotential().energy();
         }
         simulation().externalPotential().set("externalParameter",ext);
         #endif
         if (simulation().domain().isMaster()) {
            accumulator_->sample(lrf);
         }
      }
   }

   /*
   * Output results to file after simulation is completed.
   */
   void ExternalLRFAverage::output()
   {
      if (simulation().domain().isMaster()) {
         simulation().fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
         ParamComposite::writeParam(outputFile_);
         outputFile_.close();

         simulation().fileMaster().openOutputFile(outputFileName(".ave"), outputFile_);
         accumulator_->output(outputFile_);
         outputFile_.close();
      }
   }

}
