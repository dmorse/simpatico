#ifndef MCMD_MD_POTENTIAL_ENERGY_AVERAGE_CPP
#define MCMD_MD_POTENTIAL_ENERGY_AVERAGE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "MdPotentialEnergyAverage.h"        
#include <util/misc/FileMaster.h>  
#include <util/archives/Serializable_includes.h>

#include <cstdio> 

namespace McMd
{

   using namespace Util;

   /* 
   * Constructor.
   */
   MdPotentialEnergyAverage::MdPotentialEnergyAverage(MdSystem& system)
    : SystemDiagnostic<MdSystem>(system),
      outputFile_(),
      accumulator_(),
      nSamplePerBlock_(1)
   { setClassName("MdPotentialEnergyAverage"); }

   /*
   * Read parameters and initialize.
   */
   void MdPotentialEnergyAverage::readParameters(std::istream& in)
   {
      readInterval(in);
      readOutputFileName(in);
      read<int>(in,"nSamplePerBlock", nSamplePerBlock_);
      accumulator_.setNSamplePerBlock(nSamplePerBlock_);

      // If nSamplePerBlock != 0, open an output file for block averages.
      if (accumulator_.nSamplePerBlock()) {
         fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      }
      isInitialized_ = true;
   }

   /*
   * Load state from an archive.
   */
   void MdPotentialEnergyAverage::loadParameters(Serializable::IArchive& ar)
   { 
      loadInterval(ar);
      loadOutputFileName(ar);
      loadParameter<int>(ar,"nSamplePerBlock", nSamplePerBlock_);
      ar & accumulator_;

      // If nSamplePerBlock != 0, open an output file for block averages.
      if (accumulator_.nSamplePerBlock()) {
         fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      }
      isInitialized_ = true;
   }

   /*
   * Save state to archive.
   */
   void MdPotentialEnergyAverage::save(Serializable::OArchive& ar)
   { ar & *this; }

   
   /*
   * Clear accumulator.
   */
   void MdPotentialEnergyAverage::setup() 
   {  accumulator_.clear(); }
 
   /* 
   * Evaluate energy, and add to accumulator.
   */
   void MdPotentialEnergyAverage::sample(long iStep) 
   {
      if (isAtInterval(iStep))  {
         accumulator_.sample(system().potentialEnergy(), outputFile_);
      }
   }

   /*
   * Output results to file after simulation is completed.
   */
   void MdPotentialEnergyAverage::output() 
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
