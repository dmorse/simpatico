/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "StressAutoCorr.h"
//#include <util/misc/FileMaster.h>
#include <util/accumulators/AutoCorrelation.tpp>
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
      maxStageId_(10),
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
      if (simulation().domain().isMaster()) {
         // Allocate memory
         accumulator_.setParam(capacity_, maxStageId_);
      }

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
         ar >> accumulator_;
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
      if (simulation().domain().isMaster()) {
         ar << accumulator_;
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

      if (simulation().domain().isMaster()) {
         accumulator_.clear();  
      }
   }


   /*
 *    * Set actual number of molecules and clear accumulator.
 *       */
   void StressAutoCorr::setup()
   {
      if (!isInitialized_) {
         UTIL_THROW("Object not initialized.");
      }
      accumulator_.clear();
   }

   /*
   * Sample the stress tensor.
   */
   void StressAutoCorr::sample(long iStep) 
   {  
      if (isAtInterval(iStep))  {
         Simulation& sys = simulation();
         sys.computeVirialStress();
         sys.computeKineticStress();

         if (sys.domain().isMaster()) 
         {
            Tensor virial  = sys.virialStress();
            Tensor kinetic = sys.kineticStress();
            Tensor total = total.add(virial, kinetic);

            // Remove trace
            double pressure = 0.0;
            int i, j;
            for (i = 0; i < Dimension; ++i) {
               pressure += total(i,i);
            }
            pressure = pressure/double(Dimension);
            for (i = 0; i < Dimension; ++i) {
               total(i,i) -= pressure;
            }
   
            double factor = sqrt(sys.boundary().volume()/10.0);
            for (i = 0; i < Dimension; ++i) {
               for (j = 0; j < Dimension; ++j) {
                  total(i,j) *= factor;
               }
            }
   
            accumulator_.sample(total);
         }
      }
   }

   /*
   * Dump StressTensor Measurment results.
   */
   void StressAutoCorr::output() 
   {

      if (simulation().domain().isMaster()) {
         simulation().fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
         writeParam(outputFile_);
         outputFile_ << std::endl;
         outputFile_ << "bufferCapacity  " << accumulator_.bufferCapacity() << std::endl;
         outputFile_ << "nSample         " << accumulator_.nSample() << std::endl;
         outputFile_ << std::endl;
         outputFile_ << "average   " << accumulator_.average() << std::endl;
         outputFile_ << std::endl;
         outputFile_ << "Format of *.dat file" << std::endl;
         outputFile_ << "[int time in samples]  [double autocorrelation function]"
                     << std::endl;
         outputFile_ << std::endl;
         outputFile_.close();

         // Write xy autocorrelation function to data file
         simulation().fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
         accumulator_.output(outputFile_);
         outputFile_.close();
      }     
   }

}
