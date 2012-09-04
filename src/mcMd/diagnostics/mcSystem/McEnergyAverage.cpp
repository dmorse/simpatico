#ifndef MCMD_MC_ENERGY_AVERAGE_CPP
#define MCMD_MC_ENERGY_AVERAGE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "McEnergyAverage.h"                        // class header
#include <mcMd/util/FileMaster.h>  
#include <util/archives/Serializable_includes.h>


#include <cstdio> 

namespace McMd
{

   using namespace Util;

   /* 
   * Constructor.
   */
   McEnergyAverage::McEnergyAverage(McSystem& system)
    : SystemDiagnostic<McSystem>(system),
      outputFile_(),
      accumulator_(),
      nSamplePerBlock_(1),
      isInitialized_(false)
   {  setClassName("McEnergyAverage"); }

   /*
   * Read parameters and initialize.
   */
   void McEnergyAverage::readParam(std::istream& in)
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
   * Clear accumulator.
   */
   void McEnergyAverage::setup() 
   {  accumulator_.clear(); }
 
   /* 
   * Evaluate energy, and add to accumulator.
   */
   void McEnergyAverage::sample(long iStep) 
   {
      if (isAtInterval(iStep))  {
         accumulator_.sample(system().potentialEnergy(), outputFile_);
      }
   }

   /*
   * Output results to file after simulation is completed.
   */
   void McEnergyAverage::output() 
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
   
   /*
   * Save state to binary file archive.
   */
   void McEnergyAverage::save(Serializable::OArchiveType& ar)
   { ar & *this; }

   /*
   * Load state from a binary file archive.
   */
   void McEnergyAverage::load(Serializable::IArchiveType& ar)
   { ar & *this; }

}
#endif 
