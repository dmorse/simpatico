#ifndef DDMD_STRESS_TENSOR_AUTO_CORRELATION_H
#define DDMD_STRESS_TENSOR_AUTO_CORRELATION_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "StressAutoCorr.h"
//#include <util/misc/FileMaster.h>
#include <util/misc/ioUtil.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

#include <sstream>

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   StressAutoCorr::StressAutoCorr(Simulation& simulation) 
    : Analyzer(simulation),
      accumulator_(),
      capacity_(-1),
      nSample_(0),
      isInitialized_(false)
   {  setClassName("StressAutoCorr"); }

   /*
   * Read interval and outputFileName. 
   */
   void StressAutoCorr::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);
      read<int>(in,"capacity", capacity_);
      // Allocate memory
      int speciesCapacity = speciesPtr_->capacity();
      accumulator_.setParam(speciesCapacity, capacity_);

      std::string filename;
      filename  = outputFileName();
      simulation().fileMaster().openOutputFile(outputFileName(), outputFile_);

      isInitialized_ = true;
   }


   /*
   * Load internal state from an archive.
   */
   void StressAutoCorr::loadParameters(Serializable::IArchive &ar)
   {
      loadInterval(ar);
      loadOutputFileName(ar);

      MpiLoader<Serializable::IArchive> loader(*this, ar);
      loader.load(capacity_);

      if (simulation().domain().isMaster()) {
         accumulator_.loadParameters(ar);

         std::string filename;
         filename  = outputFileName();
         simulation().fileMaster().openOutputFile(outputFileName(), outputFile_);
      }
      isInitialized_ = true;
   }

   /*
   * Save internal state to an archive.
   */
   void StressAutoCorr::save(Serializable::OArchive &ar)
   {
      saveInterval(ar);
      saveOutputFileName(ar);
      ar << nSample_;

      if (simulation().domain().isMaster()) {
         ar << Accumulator_;
      }

   }

  
   /*
   * Read interval and outputFileName. 
   */
   void StressAutoCorr::clear() 
   {
      if (!isInitialized_) {
         UTIL_THROW("Error: Object not initialized");
      }
      nSample_ = 0;
      if (simulation().domain().isMaster()) {
         Accumulator_.clear();  
      }
   }

   /*
   * Sample the stress tensor.
   */
   void StressAutoCorr::sample(long iStep) 
   {  
      double pressure;
      double strAutoCorr; 
      if (isAtInterval(iStep))  {
         Simulation& sys = simulation();
         sys.computeVirialStress();
         sys.computeKineticStress();
         if (sys.domain().isMaster()) {
            Tensor virial  = sys.virialStress();
            Tensor kinetic = sys.kineticStress();
            Tensor total = total.add(virial, kinetic);
            pressure = sys.kineticPressure()+sys.virialPressure();
            strAutoCorr = (virial(0,0)-pressure)
             
            accumulator_.sample(virial(0,0));
         }
         ++nSample_;
      }
   }

   /*
   * Dump StressTensor Measurment results.
   */
   void StressAutoCorr::output() 
   {

      if (simulation().domain().isMaster()) {
        
         outputFile_ << 
                     "Sxx=" << Dbl(sxxAccumulator_.average(), 17)<< "  +-  " << Dbl(sxxAccumulator_.error(), 9, 8) << "\n" << 
                     "Sxy=" << Dbl(sxyAccumulator_.average(), 17)<< "  +-  " << Dbl(sxyAccumulator_.error(), 9, 8) << "\n" << 
                     "Sxz=" << Dbl(sxzAccumulator_.average(), 17)<< "  +-  " << Dbl(sxzAccumulator_.error(), 9, 8) << "\n" << 
                     "Syx=" << Dbl(syxAccumulator_.average(), 17)<< "  +-  " << Dbl(syxAccumulator_.error(), 9, 8) << "\n" << 
                     "Syy=" << Dbl(syyAccumulator_.average(), 17)<< "  +-  " << Dbl(syyAccumulator_.error(), 9, 8) << "\n" << 
                     "Syz=" << Dbl(syzAccumulator_.average(), 17)<< "  +-  " << Dbl(syzAccumulator_.error(), 9, 8) << "\n" << 
                     "Szx=" << Dbl(szxAccumulator_.average(), 17)<< "  +-  " << Dbl(szxAccumulator_.error(), 9, 8) << "\n" << 
                     "Szy=" << Dbl(szyAccumulator_.average(), 17)<< "  +-  " << Dbl(szyAccumulator_.error(), 9, 8) << "\n" << 
                     "Szz=" << Dbl(szzAccumulator_.average(), 17)<< "  +-  " << Dbl(szzAccumulator_.error(), 9, 8) << "\n" << 
                     std::endl;
         }
   }

}
#endif  
