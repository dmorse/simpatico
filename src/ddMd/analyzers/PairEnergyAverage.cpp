#ifndef DDMD_PAIR_ENERGY_AVERAGE_CPP
#define DDMD_PAIR_ENERGY_AVERAGE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "PairEnergyAverage.h"
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
   PairEnergyAverage::PairEnergyAverage(Simulation& simulation) 
    : Analyzer(simulation),
      outputFile_(),
      pairs_(),
      accumulator_(),
      nSamplePerBlock_(1),
      isInitialized_(false)
   {  setClassName("PairEnergyAverage"); }

   /*
   * Read interval and outputFileName. 
   */
   void PairEnergyAverage::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);

      pairs_.allocate(2);
      readDArray<int>(in, "pairs", pairs_, 2);

      read<int>(in,"nSamplePerBlock", nSamplePerBlock_);
      accumulator_.setNSamplePerBlock(nSamplePerBlock_);
      accumulator_.clear();

      if (sys.domain().isMaster()) {
      simulation().fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      }

      isInitialized_ = true;
   }

   /*
   * Load internal state from an archive.
   */
   void PairEnergyAverage::loadParameters(Serializable::IArchive &ar)
   {
      loadInterval(ar);
      MpiLoader<Serializable::IArchive> loader(*this, ar);
      loader.load(nSamplePerBlock_);
      ar & accumulator_;

      if (nSamplePerBlock_ != accumulator_.nSamplePerBlock()) {
         UTIL_THROW("Inconsistent values of nSamplePerBlock");
      }

      isInitialized_ = true;
   }

   /*
   * Save internal state to an archive.
   */
   void PairEnergyAverage::save(Serializable::OArchive &ar)
   { ar & *this; }

   /*
   * Reset nSample.
   */
   void PairEnergyAverage::clear() 
   {  accumulator_.clear();  }

   /*
   * Dump configuration to file
   */
   void PairEnergyAverage::sample(long iStep) 
   {
      if (isAtInterval(iStep))  {
         Simulation& sys = simulation();
         sys.computePairEnergies();
         if (sys.domain().isMaster()) {
            DMatrix<double> pair = sys.pairEnergies();
            for (int i = 0; i < simulation().nAtomType(); ++i){
               for (int j = 0; j < simulation().nAtomType(); ++j){
                  pair(i,j) = 0.5*( pair(i,j)+pair(j,i) );
                  pair(j,i) = pair(i,j);
               }
            }

            outputFile_ << Dbl(pair(pairs_[0],pairs_[1]), 15);
            accumulator_.sample(pair(pairs_[0],pairs_[1]));
         }
      }
   }

   /*
   * Output results to file after simulation is completed.
   */
   void PairEnergyAverage::output()
   {
      Simulation& sys = simulation();
      if (sys.domain().isMaster()) {
      outputFile_.close();
      simulation().fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
      ParamComposite::writeParam(outputFile_);
      outputFile_.close();
      
      simulation().fileMaster().openOutputFile(outputFileName(".ave"), outputFile_);
      accumulator_.output(outputFile_);
      outputFile_.close();
      }
   }

}
#endif
