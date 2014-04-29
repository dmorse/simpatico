#ifndef DDMD_EXTERNAL_POTENTIAL_AVERAGE_CPP
#define DDMD_EXTERNAL_POTENTIAL_AVERAGE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "ExternalPotentialAverage.h"
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
   ExternalPotentialAverage::ExternalPotentialAverage(Simulation& simulation) 
    : Analyzer(simulation),
      outputFile_(),
      accumulator_(),
      nSamplePerBlock_(1),
      isInitialized_(false)
   {  setClassName("ExternalPotentialAverage"); }

   /*
   * Read interval and outputFileName. 
   */
   void ExternalPotentialAverage::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);
      read<int>(in,"nSamplePerBlock", nSamplePerBlock_);

      accumulator_.setNSamplePerBlock(nSamplePerBlock_);

      isInitialized_ = true;
   }

   /*
   * Load internal state from an archive.
   */
   void ExternalPotentialAverage::loadParameters(Serializable::IArchive &ar)
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
   void ExternalPotentialAverage::save(Serializable::OArchive &ar)
   {
      ar & *this;
   }

   /*
   * Reset nSample.
   */
   void ExternalPotentialAverage::clear() 
   {  accumulator_.clear();  }

   /*
   * Dump configuration to file
   */
   void ExternalPotentialAverage::sample(long iStep) 
   {
      if (isAtInterval(iStep))  {
         Simulation& sys = simulation();
         if (sys.domain().isMaster()) {
            double potential = 0.0;

            #ifdef INTER_EXTERNAL
            if (sys.nExternalType()) {
               double external  = sys.externalPotential().energy();
               potential += external;
            }
            #endif

            accumulator_.sample(potential);
         }
      }
   }

   /*
   * Output results to file after simulation is completed.
   */
   void ExternalPotentialAverage::output()
   {
      Simulation& sys = simulation();
      if (sys.domain().isMaster()) {
      // Write parameters to file
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
