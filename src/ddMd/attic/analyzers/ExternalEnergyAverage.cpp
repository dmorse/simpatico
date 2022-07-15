/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ExternalEnergyAverage.h"
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
   ExternalEnergyAverage::ExternalEnergyAverage(Simulation& simulation) 
    : Analyzer(simulation),
      outputFile_(),
      accumulator_(NULL),
      nSamplePerBlock_(1),
      isInitialized_(false)
   {  setClassName("ExternalEnergyAverage"); }

   /*
   * Destructor.
   */
   ExternalEnergyAverage::~ExternalEnergyAverage() 
   {  
      if(accumulator_ != NULL) {
         delete accumulator_;
      }
   }

   /*
   * Read interval and outputFileName. 
   */
   void ExternalEnergyAverage::readParameters(std::istream& in) 
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
   void ExternalEnergyAverage::loadParameters(Serializable::IArchive &ar)
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
   void ExternalEnergyAverage::save(Serializable::OArchive &ar)
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
   void ExternalEnergyAverage::clear() 
   {   
      if (simulation().domain().isMaster()){ 
         accumulator_->clear();
      }
   }

   /*
   * Dump configuration to file
   */
   void ExternalEnergyAverage::sample(long iStep) 
   {
      if (isAtInterval(iStep))  {
         double external = 0.0;
         simulation().computePotentialEnergies();
         #ifdef SIMP_EXTERNAL
         if (simulation().hasExternal()) {
            external =  simulation().externalPotential().energy();
         }
         #endif
         if (simulation().domain().isMaster()) {
            accumulator_->sample(external);
         }
      }
   }

   /*
   * Output results to file after simulation is completed.
   */
   void ExternalEnergyAverage::output()
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
