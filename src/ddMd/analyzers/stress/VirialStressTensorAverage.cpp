/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "VirialStressTensorAverage.h"
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
   VirialStressTensorAverage::VirialStressTensorAverage(Simulation& simulation) 
    : Analyzer(simulation),
      sxxAccumulator_(),
      sxyAccumulator_(),
      sxzAccumulator_(),
      syxAccumulator_(),
      syyAccumulator_(),
      syzAccumulator_(),
      szxAccumulator_(),
      szyAccumulator_(),
      szzAccumulator_(),
      nSamplePerBlock_(1),
      isInitialized_(false)
   {  setClassName("VirialStressTensorAverage"); }

   /*
   * Read interval and outputFileName. 
   */
   void VirialStressTensorAverage::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);
      read<int>(in,"nSamplePerBlock", nSamplePerBlock_);
      sxxAccumulator_.setNSamplePerBlock(nSamplePerBlock_);
      sxyAccumulator_.setNSamplePerBlock(nSamplePerBlock_);
      sxzAccumulator_.setNSamplePerBlock(nSamplePerBlock_);
      syxAccumulator_.setNSamplePerBlock(nSamplePerBlock_);
      syyAccumulator_.setNSamplePerBlock(nSamplePerBlock_);
      syzAccumulator_.setNSamplePerBlock(nSamplePerBlock_);
      szxAccumulator_.setNSamplePerBlock(nSamplePerBlock_);
      szyAccumulator_.setNSamplePerBlock(nSamplePerBlock_);
      szzAccumulator_.setNSamplePerBlock(nSamplePerBlock_);
      simulation().fileMaster().openOutputFile(outputFileName(), outputFile_);
      isInitialized_ = true;
   }


   /*
   * Load internal state from an archive.
   */
   void VirialStressTensorAverage::loadParameters(Serializable::IArchive &ar)
   {
      loadInterval(ar);
      loadOutputFileName(ar);
      loadParameter<int>(ar,"nSamplePerBlock", nSamplePerBlock_);
      if (simulation().domain().isMaster()) {
         #if 0
         sxxAccumulator_.loadParameters(ar);
         sxyAccumulator_.loadParameters(ar);
         sxzAccumulator_.loadParameters(ar);
         syxAccumulator_.loadParameters(ar);
         syyAccumulator_.loadParameters(ar);
         syzAccumulator_.loadParameters(ar);
         szxAccumulator_.loadParameters(ar);
         szyAccumulator_.loadParameters(ar);
         szzAccumulator_.loadParameters(ar);
         #endif
         ar >> sxxAccumulator_;
         ar >> sxyAccumulator_;
         ar >> sxzAccumulator_;
         ar >> syxAccumulator_;
         ar >> syyAccumulator_;
         ar >> syzAccumulator_;
         ar >> szxAccumulator_;
         ar >> szyAccumulator_;
         ar >> szzAccumulator_;
         simulation().fileMaster().openOutputFile(outputFileName(), outputFile_);
      }
      isInitialized_ = true;
   }

   /*
   * Save internal state to an archive.
   */
   void VirialStressTensorAverage::save(Serializable::OArchive &ar)
   {
      assert(simulation().domain().isMaster());
      saveInterval(ar);
      saveOutputFileName(ar);
      ar << nSamplePerBlock_;
      ar << sxxAccumulator_;
      ar << sxyAccumulator_;
      ar << sxzAccumulator_;
      ar << syxAccumulator_;
      ar << syyAccumulator_;
      ar << syzAccumulator_;
      ar << szxAccumulator_;
      ar << szyAccumulator_;
      ar << szzAccumulator_;
   }

  
   /*
   * Read interval and outputFileName. 
   */
   void VirialStressTensorAverage::clear() 
   {
      if (!isInitialized_) {
         UTIL_THROW("Error: Object not initialized");
      }
      nSample_ = 0;
      if (simulation().domain().isMaster()) {
         sxxAccumulator_.clear();
         sxyAccumulator_.clear();
         sxzAccumulator_.clear();
         syxAccumulator_.clear();
         syyAccumulator_.clear();
         syzAccumulator_.clear();
         szxAccumulator_.clear();
         szyAccumulator_.clear();
         szzAccumulator_.clear();       
      }
   }

   /*
   * Sample the stress tensor.
   */
   void VirialStressTensorAverage::sample(long iStep) 
   {
      if (isAtInterval(iStep))  {
         Simulation& sys = simulation();
         sys.computeVirialStress();
         sys.computeKineticStress();
         if (sys.domain().isMaster()) {
            Tensor virial  = sys.virialStress();
            Tensor kinetic = sys.kineticStress();
            Tensor total = total.add(virial, kinetic); 
            sxxAccumulator_.sample(virial(0,0));
            sxyAccumulator_.sample(virial(0,1));
            sxzAccumulator_.sample(virial(0,2));
            syxAccumulator_.sample(virial(1,0));
            syyAccumulator_.sample(virial(1,1));
            syzAccumulator_.sample(virial(1,2));
            szxAccumulator_.sample(virial(2,0));
            szyAccumulator_.sample(virial(2,1));
            szzAccumulator_.sample(virial(2,2));
         }
         ++nSample_;
      }
   }

   /*
   * Dump StressTensor Measurment results.
   */
   void VirialStressTensorAverage::output() 
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
