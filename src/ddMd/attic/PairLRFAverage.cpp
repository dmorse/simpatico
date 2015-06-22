/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "PairLRFAverage.h"
#include <ddMd/potentials/pair/PairPotential.h>
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
   PairLRFAverage::PairLRFAverage(Simulation& simulation) 
    : Analyzer(simulation),
      outputFile_(),
      pairs_(),
      accumulator_(NULL),
      nSamplePerBlock_(1),
      isInitialized_(false)
   {  setClassName("PairLRFAverage"); }

   /*
   * Destructor.
   */
   PairLRFAverage::~PairLRFAverage() 
   {
      if(simulation().domain().isMaster()) {
         delete accumulator_;
      }   
   }

   /*
   * Read interval and outputFileName. 
   */
   void PairLRFAverage::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);
      pairs_.allocate(2);
      readDArray<int>(in, "pairs", pairs_, 2);
      read<int>(in,"nSamplePerBlock", nSamplePerBlock_);
      if(simulation().domain().isMaster()) {
         accumulator_ = new Average;
         accumulator_->setNSamplePerBlock(nSamplePerBlock_);
      }

      isInitialized_ = true;
   }

   /*
   * Load internal state from an archive.
   */
   void PairLRFAverage::loadParameters(Serializable::IArchive &ar)
   {
      loadInterval(ar);
      loadOutputFileName(ar);
      pairs_.allocate(2);
      loadDArray<int>(ar, "pairs", pairs_, 2);
      loadParameter<int>(ar,"nSamplePerBlock", nSamplePerBlock_);
      if (simulation().domain().isMaster()) {
         accumulator_ = new Average;
         ar >> *accumulator_;
         if (nSamplePerBlock_ != accumulator_->nSamplePerBlock()) {
            UTIL_THROW("Inconsistent values of nSamplePerBlock");
         }
      }
      isInitialized_ = true;
   }

   /*
   * Save internal state to an archive.
   */
   void PairLRFAverage::save(Serializable::OArchive &ar)
   {       
      saveInterval(ar);
      saveOutputFileName(ar);
      ar << pairs_;
      ar << nSamplePerBlock_;
      ar << *accumulator_;
   }

   /*
   * Reset nSample.
   */
   void PairLRFAverage::clear() 
   {  
      if (simulation().domain().isMaster()){
         accumulator_->clear();
      } 
   }

   /*
   * Dump configuration to file
   */
   void PairLRFAverage::sample(long iStep) 
   {
      if (isAtInterval(iStep))  {
         double epsilon11 = simulation().pairPotential().get("epsilon", 1, 1);
         double epsilon12 = simulation().pairPotential().get("epsilon", 1, 2);
         double epsilon21 = simulation().pairPotential().get("epsilon", 2, 1);
         double epsilon22 = simulation().pairPotential().get("epsilon", 2, 2);
         simulation().pairPotential().set("epsilon", 1, 1, 1.0);
         simulation().pairPotential().set("epsilon", 1, 2, 1.0);
         simulation().pairPotential().set("epsilon", 2, 1, 1.0);
         simulation().pairPotential().set("epsilon", 2, 2, 1.0);
         simulation().computePairEnergies();
         if (simulation().domain().isMaster()) {
            DMatrix<double> pair = simulation().pairEnergies();
            for (int i = 0; i < simulation().nAtomType(); ++i){
               for (int j = 0; j < simulation().nAtomType(); ++j){
                  pair(i,j) = 0.5*( pair(i,j)+pair(j,i) );
                  pair(j,i) = pair(i,j);
               }
            }
            accumulator_->sample(pair(pairs_[0],pairs_[1]));
         }
         simulation().pairPotential().set("epsilon", 1, 1, epsilon11);
         simulation().pairPotential().set("epsilon", 1, 2, epsilon12);
         simulation().pairPotential().set("epsilon", 2, 1, epsilon21);
         simulation().pairPotential().set("epsilon", 2, 2, epsilon22);
      }
   }

   /*
   * Output results to file after simulation is completed.
   */
   void PairLRFAverage::output()
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
