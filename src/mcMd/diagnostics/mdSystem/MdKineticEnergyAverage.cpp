#ifndef MCMD_MD_KINETIC_ENERGY_AVERAGE_CPP
#define MCMD_MD_KINETIC_ENERGY_AVERAGE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "MdKineticEnergyAverage.h"        // class header
#include <util/misc/FileMaster.h>  
#include <util/archives/Serializable_includes.h>

#include <cstdio> 

namespace McMd
{

   using namespace Util;

   /* 
   * Constructor.
   */
   MdKineticEnergyAverage::MdKineticEnergyAverage(MdSystem& system)
    : SystemDiagnostic<MdSystem>(system),
      outputFile_(),
      accumulator_(),
      nSamplePerBlock_(1)
   { setClassName("MdKineticEnergyAverage"); }

   /*
   * Read parameters and initialize.
   */
   void MdKineticEnergyAverage::readParameters(std::istream& in)
   {
      readInterval(in);
      readOutputFileName(in);
      read<int>(in, "nSamplePerBlock", nSamplePerBlock_);

      accumulator_.setNSamplePerBlock(nSamplePerBlock_);

      // If nSamplePerBlock != 0, open an output file for block averages.
      if (accumulator_.nSamplePerBlock()) {
         fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      }
      isInitialized_ = true;
   }

   /*
   * Load internal state from archive.
   */
   void MdKineticEnergyAverage::loadParameters(Serializable::IArchive &ar)
   {
      loadInterval(ar);
      loadOutputFileName(ar);
      loadParameter<int>(ar, "nSamplePerBlock", nSamplePerBlock_);
      ar & accumulator_;

      // If nSamplePerBlock != 0, open an output file for block averages.
      if (accumulator_.nSamplePerBlock()) {
         fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      }
      isInitialized_ = true;
   }

   /*
   * Save internal state to archive.
   */
   void MdKineticEnergyAverage::save(Serializable::OArchive &ar)
   { ar & *this; }
   
   /*
   * Clear accumulator.
   */
   void MdKineticEnergyAverage::setup() 
   {  accumulator_.clear(); }
 
   /* 
   * Evaluate energy, and add to accumulator.
   */
   void MdKineticEnergyAverage::sample(long iStep) 
   {
      if (isAtInterval(iStep))  {
         accumulator_.sample(system().kineticEnergy(), outputFile_);
      }
   }

   /*
   * Output results to file after simulation is completed.
   */
   void MdKineticEnergyAverage::output() 
   { 
      // If outputFile_ was used to write block averages, close it.
      if (accumulator_.nSamplePerBlock()) {
         outputFile_.close();
      }

      fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
      writeParam(outputFile_); 
      outputFile_.close();

      fileMaster().openOutputFile(outputFileName(".ave"), outputFile_);
      accumulator_.output(outputFile_); 
      outputFile_.close();

   }

}
#endif 
