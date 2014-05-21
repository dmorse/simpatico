#ifndef DDMD_STRESS_TENSOR_AUTO_CORRELATION_CPP
#define DDMD_STRESS_TENSOR_AUTO_CORRELATION_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "StressAutoCorrelation.h"
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
   StressAutoCorrelation::StressAutoCorrelation(Simulation& simulation) 
    : Analyzer(simulation),
      accumulator_(),
      capacity_(-1),
      isInitialized_(false)
   {  setClassName("StressAutoCorrelation"); }

   /*
   * Read interval and outputFileName. 
   */
   void StressAutoCorrelation::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);
      read<int>(in,"capacity", capacity_);
      // Allocate memory
      // Autocorrelation accumulator allocation is so to work for
      // a symmetric matrix, 6 is number of independent elements in 
      // a symmetric matrix.
      accumulator_.setParam(6, capacity_);

      isInitialized_ = true;
   }


   /*
   * Load internal state from an archive.
   */
   void StressAutoCorrelation::loadParameters(Serializable::IArchive &ar)
   {
      loadInterval(ar);
      loadOutputFileName(ar);

      MpiLoader<Serializable::IArchive> loader(*this, ar);
      loader.load(capacity_);

      if (simulation().domain().isMaster()) {
         accumulator_.loadParameters(ar);
      }

      isInitialized_ = true;
   }

   /*
   * Save internal state to an archive.
   */
   void StressAutoCorrelation::save(Serializable::OArchive &ar)
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
   void StressAutoCorrelation::clear() 
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
   void StressAutoCorrelation::setup()
   {
      if (!isInitialized_) {
         UTIL_THROW("Object not initialized.");
      }

       // Set number of molecules and clear accumulator
       accumulator_.setNEnsemble(6);
       accumulator_.clear();
   }

   /*
   * Sample the stress tensor.
   */
   void StressAutoCorrelation::sample(long iStep) 
   {  
      DArray<double> elements;
      elements.allocate(9);
      double pressure;
      double temprature;
 
      if (isAtInterval(iStep))  {
         Simulation& sys = simulation();
         sys.computeVirialStress();
         sys.computeKineticStress();
         sys.computeKineticEnergy();
         simulation().atomStorage().computeNAtomTotal(simulation().domain().communicator());

         if (sys.domain().isMaster()) {
            Tensor virial  = sys.virialStress();
            Tensor kinetic = sys.kineticStress();
            Tensor total = total.add(virial, kinetic);
            pressure = sys.kineticPressure()+sys.virialPressure();
            double ndof = simulation().atomStorage().nAtomTotal()*3;
            //temprature = sys.kineticEnergy()*2.0/ndof;
            temprature = 1;

             
            elements[0] = (total(0,0) - pressure / 3.0) / (10.0 * temprature);
            elements[4] = (total(1,1) - pressure / 3.0) / (10.0 * temprature);
            elements[8] = (total(2,2) - pressure / 3.0) / (10.0 * temprature);
            elements[1] = (total(0,1) + total(1,0)) / 2.0 / (10.0 * temprature);
            elements[2] = (total(0,2) + total(2,0)) / 2.0 / (10.0 * temprature);
            elements[5] = (total(1,2) + total(2,1)) / 2.0 / (10.0 * temprature);
            elements[3] = elements[1];
            elements[6] = elements[2];
            elements[7] = elements[5];
             
          
            accumulator_.sample(elements);
         }
      }
   }

   /*
   * Dump StressTensor Measurment results.
   */
   void StressAutoCorrelation::output() 
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
         simulation().fileMaster().openOutputFile(outputFileName(".corr"), outputFile_);
         accumulator_.output(outputFile_);
         outputFile_.close();
      }     
   }

}
#endif  
