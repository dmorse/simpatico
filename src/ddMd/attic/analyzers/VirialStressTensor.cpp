/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "VirialStressTensor.h"
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
   VirialStressTensor::VirialStressTensor(Simulation& simulation) 
    : Analyzer(simulation),
      nSample_(0),
      isInitialized_(false)
   {  setClassName("VirialStressTensor"); }

   /*
   * Read interval and outputFileName. 
   */
   void VirialStressTensor::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);
     
      if(simulation().domain().isMaster()) {
         std::string filename;
         filename  = outputFileName();
         simulation().fileMaster().openOutputFile(outputFileName(), outputFile_);
      }

      isInitialized_ = true;
   }


   /*
   * Load internal state from an archive.
   */
   void VirialStressTensor::loadParameters(Serializable::IArchive &ar)
   {
      loadInterval(ar);
      loadOutputFileName(ar);

      MpiLoader<Serializable::IArchive> loader(*this, ar);
      loader.load(nSample_);

      if (simulation().domain().isMaster()) {
         std::string filename;
         filename  = outputFileName();
         simulation().fileMaster().openOutputFile(outputFileName(), outputFile_);
      }

      isInitialized_ = true;
   }

   /*
   * Save internal state to an archive.
   */
   void VirialStressTensor::save(Serializable::OArchive &ar)
   {
      saveInterval(ar);
      saveOutputFileName(ar);
      ar << nSample_;
   }

  
   /*
   * Read interval and outputFileName. 
   */
   void VirialStressTensor::clear() 
   {
      if (!isInitialized_) {
         UTIL_THROW("Error: Object not initialized");
      }
      nSample_ = 0;
   }

   /*
   * Sample the stress tensor.
   */
   void VirialStressTensor::sample(long iStep) 
   {
      if (isAtInterval(iStep))  {
         simulation().computeVirialStress();
         simulation().computeKineticStress();
         if (simulation().domain().isMaster()) {
            Tensor virial  = simulation().virialStress();
            Tensor kinetic = simulation().kineticStress();
            Tensor total;
            total.add(virial, kinetic); 
            outputFile_ << virial(0,0) << "\t" << virial(1,1) << "\t" << virial(2,2)
            << total(0,0) << "\t" << total(1,1) << "\t" << total(2,2) << "\n"; 
         }
         
         ++nSample_;
      }
   }

}
