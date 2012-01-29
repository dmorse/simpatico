#ifndef MD_KINETIC_ENERGY_AVERAGE_CPP
#define MD_KINETIC_ENERGY_AVERAGE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "MdKineticEnergyAverage.h"        // class header
#include <mcMd/util/FileMaster.h>  
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
   {}

   /*
   * Read parameters and initialize.
   */
   void MdKineticEnergyAverage::readParam(std::istream& in)
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
   void MdKineticEnergyAverage::initialize() 
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

   /*
   * Save state to binary file archive.
   */
   void MdKineticEnergyAverage::save(Serializable::OArchiveType& ar)
   { ar & *this; }

   /*
   * Load state from a binary file archive.
   */
   void MdKineticEnergyAverage::load(Serializable::IArchiveType& ar)
   { ar & *this; }

   
}
#endif 
