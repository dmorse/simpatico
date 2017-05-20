/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "OutputStressTensor.h"
#include <util/misc/FileMaster.h>
#include <util/mpi/MpiLoader.h>
#include <util/space/Tensor.h>
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
   OutputStressTensor::OutputStressTensor(Simulation& simulation) 
    : Analyzer(simulation),
      nSample_(0),
      isInitialized_(false)
   {  setClassName("OutputStressTensor"); }

   /*
   * Read interval and outputFileName. 
   */
   void OutputStressTensor::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);

      #if 0
      std::string filename;
      filename  = outputFileName();
      simulation().fileMaster().openOutputFile(filename, outputFile_);
      #endif
      isInitialized_ = true;
   }

   /*
   * Load internal state from an archive.
   */
   void OutputStressTensor::loadParameters(Serializable::IArchive &ar)
   {
      loadInterval(ar);
      loadOutputFileName(ar);

      MpiLoader<Serializable::IArchive> loader(*this, ar);
      loader.load(nSample_);

      #if 0
      std::string filename;
      filename  = outputFileName();
      simulation().fileMaster().openOutputFile(filename, outputFile_);
      #endif
      isInitialized_ = true;
   }

   /*
   * Save internal state to an archive.
   */
   void OutputStressTensor::save(Serializable::OArchive &ar)
   {
      saveInterval(ar);
      saveOutputFileName(ar);
      ar << nSample_;
   }

  
   /*
   * Read interval and outputFileName. 
   */
   void OutputStressTensor::clear() 
   {  nSample_ = 0; }

   /*
   * Open outputfile
   */ 
   void OutputStressTensor::setup()
   {
      if (simulation().domain().isMaster()) {
         std::string filename;
         filename  = outputFileName();
         simulation().fileMaster().openOutputFile(filename, outputFile_);
      }
   }

   /*
   * Dump configuration to file
   */
   void OutputStressTensor::sample(long iStep) 
   {
      if (isAtInterval(iStep))  {
         Simulation& sys = simulation();
         sys.computeVirialStress();
         sys.computeKineticStress();
         if (sys.domain().isMaster()) {
            Tensor virial  = sys.virialStress();
            Tensor kinetic = sys.kineticStress();
            Tensor total;
            total.add(virial, kinetic);
            outputFile_ << Int(iStep, 10)
                        << Dbl(virial(0,0), 20)
                        << Dbl(virial(0,1), 20)
                        << Dbl(virial(0,2), 20)
                        << Dbl(virial(1,0), 20)
                        << Dbl(virial(1,1), 20)
                        << Dbl(virial(1,2), 20)
                        << Dbl(virial(2,0), 20)
                        << Dbl(virial(2,1), 20)
                        << Dbl(virial(2,2), 20)
                        << Dbl(kinetic(0,0), 20)
                        << Dbl(kinetic(0,1), 20)
                        << Dbl(kinetic(0,2), 20)
                        << Dbl(kinetic(1,0), 20)
                        << Dbl(kinetic(1,1), 20)
                        << Dbl(kinetic(1,2), 20)
                        << Dbl(kinetic(2,0), 20)
                        << Dbl(kinetic(2,1), 20)
                        << Dbl(kinetic(2,2), 20)
                        << Dbl(total(0,0), 20)
                        << Dbl(total(0,1), 20)
                        << Dbl(total(0,2), 20)
                        << Dbl(total(1,0), 20)
                        << Dbl(total(1,1), 20)
                        << Dbl(total(1,2), 20)
                        << Dbl(total(2,0), 20)
                        << Dbl(total(2,1), 20)
                        << Dbl(total(2,2), 20)
                        << std::endl;
         }

         ++nSample_;
      }
   }

}
