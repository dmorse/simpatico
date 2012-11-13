#ifndef MCMD_MC_STRESS_AUTO_CORR_CPP
#define MCMD_MC_STRESS_AUTO_CORR_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "McStressAutoCorr.h"
#include <mcMd/simulation/Simulation.h>
#include <mcMd/mcSimulation/McSystem.h>
#include <util/misc/FileMaster.h>
#include <util/archives/Serializable_includes.h>
#include <util/global.h>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   McStressAutoCorr::McStressAutoCorr(McSystem& system) 
    : SystemDiagnostic<McSystem>(system),
      outputFile_(),
      accumulator_(),
      capacity_(-1),
      isInitialized_(false)
   {  setClassName("McStressAutoCorr"); }

   /*
   * Read parameters from file, and allocate data array.
   */
   void McStressAutoCorr::readParameters(std::istream& in) 
   {

      // Read interval and parameters for AutoCorrArray
      readInterval(in);
      readOutputFileName(in);
      read<int>(in, "capacity", capacity_);

      // Validate input
      if (capacity_ <= 0) {
         UTIL_THROW("Negative capacity");
      }

      accumulator_.setParam(capacity_);
      fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      isInitialized_ = true;
   }

   /*
   * Load internal state from archive.
   */
   void McStressAutoCorr::loadParameters(Serializable::IArchive &ar)
   {
      Diagnostic::loadParameters(ar);
      loadParameter<int>(ar, "capacity", capacity_);
      ar & accumulator_;

      // Validate input
      if (capacity_ <= 0) {
         UTIL_THROW("Negative capacity");
      }
      if (capacity_ != accumulator_.bufferCapacity()) {
         UTIL_THROW("Inconsistent values for capacity");
      }

      fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      isInitialized_ = true;
   }

   /*
   * Save internal state to archive (use serialize method).
   */
   void McStressAutoCorr::save(Serializable::OArchive &ar)
   {  ar & *this; }
   

   /*
   * Clear accumulator.
   */
   void McStressAutoCorr::setup()
   {  accumulator_.clear(); }

   /*
   * Evaluate shear stress autocorrelation function.
   */
   void McStressAutoCorr::sample(long iStep) 
   { 
      if (isAtInterval(iStep))  {
         Tensor stress;
         system().computeStress(stress);
         double pressure = stress.trace()/double(Dimension);
         for (int i = 0; i < Dimension; ++i) {
            stress(i, i) -= pressure;
         }
         stress.symmetrize();
         double volume = system().boundary().volume();
         stress *= sqrt(volume);
         outputFile_ << stress << std::endl;
         accumulator_.sample(stress);
      } 
   }

   /*
   * Output results to file after simulation is completed.
   */
   void McStressAutoCorr::output() 
   {  
      outputFile_.close();

      // Echo parameters to log file
      fileMaster().openOutputFile(outputFileName(), outputFile_);
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

      // Write autocorrelation function to data file
      fileMaster().openOutputFile(outputFileName("_corr.dat"), outputFile_);
      accumulator_.output(outputFile_); 
      outputFile_.close();
   }

}

#endif 
