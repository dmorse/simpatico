/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "StressAutoCorr.h"
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
      accumulatorPtr_(0),
      bufferCapacity_(-1),
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
      read<int>(in,"bufferCapacity", bufferCapacity_);
      if (simulation().domain().isMaster()) {
         accumulatorPtr_ = new AutoCorrelation<Tensor, double>;
         accumulatorPtr_->setParam(bufferCapacity_, maxStageId_);
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
      loadParameter(ar, "bufferCapacity", bufferCapacity_);

      if (simulation().domain().isMaster()) {
         accumulatorPtr_ = new AutoCorrelation<Tensor, double>;
         ar >> *accumulatorPtr_;
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
      ar & bufferCapacity_;
      if (simulation().domain().isMaster()) {
         if (!accumulatorPtr_) {
            UTIL_THROW("Null accumulatorPtr_ on master");
         }
         ar << *accumulatorPtr_;
      }
   }
  
   /*
   * Clear accumulator.
   */
   void StressAutoCorr::clear() 
   {
      if (!isInitialized_) {
         UTIL_THROW("Error: Object not initialized");
      }
      if (simulation().domain().isMaster()) {
         if (!accumulatorPtr_) {
            UTIL_THROW("Null accumulatorPtr_ on master");
         }
         accumulatorPtr_->clear();  
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
      clear();
   }

   /*
   * Sample the stress tensor.
   */
   void StressAutoCorr::sample(long iStep) 
   {  
      if (isAtInterval(iStep))  {
         Simulation& sim = simulation();
         sim.computeVirialStress();
         sim.computeKineticStress();

         if (sim.domain().isMaster()) {
            if (!accumulatorPtr_) {
               UTIL_THROW("Null accumulatorPtr_ on master");
            }

            Tensor virial  = sim.virialStress();
            Tensor kinetic = sim.kineticStress();
            Tensor total;
            total.add(virial, kinetic);

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
   
            double factor = sqrt(sim.boundary().volume()/10.0);
            for (i = 0; i < Dimension; ++i) {
               for (j = 0; j < Dimension; ++j) {
                  total(i,j) *= factor;
               }
            }
   
            accumulatorPtr_->sample(total);
         }
      }
   }

   /*
   * Dump StressTensor Measurment results.
   */
   void StressAutoCorr::output() 
   {
      if (simulation().domain().isMaster()) {
         if (!accumulatorPtr_) {
            UTIL_THROW("Null accumulatorPtr_ on master");
         }

         simulation().fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
         writeParam(outputFile_);
         outputFile_ << std::endl;
         outputFile_ << "bufferCapacity  " << accumulatorPtr_->bufferCapacity() << std::endl;
         outputFile_ << "nSample         " << accumulatorPtr_->nSample() << std::endl;
         outputFile_ << std::endl;
         outputFile_ << "Format of *.dat file" << std::endl;
         outputFile_ << "[int time delay (samples)]  [double autocorrelation function]"
                     << std::endl;
         outputFile_ << std::endl;
         outputFile_.close();

         // Write xy autocorrelation function to data file
         simulation().fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
         accumulatorPtr_->output(outputFile_);
         outputFile_.close();
      }
   }

}
