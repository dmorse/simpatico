#ifndef VOLUME_AVERAGE_CPP
#define VOLUME_AVERAGE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "VolumeAverage.h"
#include <mcMd/simulation/Simulation.h>
#include <mcMd/species/Species.h>
#include <mcMd/boundary/Boundary.h>
#include <util/space/Dimension.h>
#include <util/archives/Serializable_includes.h>

#include <util/global.h>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   VolumeAverage::VolumeAverage(System& system) 
    : SystemDiagnostic<System>(system),
      outputFile_(),
      accumulator_(),
      nSamplePerBlock_(-1),
      isInitialized_(false)
   {}

   /*
   * Read parameters from file, and allocate data array.
   */
   void VolumeAverage::readParam(std::istream& in) 
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
   * Initialize at beginning of simulation.
   */
   void VolumeAverage::initialize() 
   {
      if (!isInitialized_) {
         UTIL_THROW("Error: object is not initialized");
      }
      accumulator_.clear();

   }

   /*
   * Evaluate end-to-end vectors of all chains, add to ensemble.
   */
   void VolumeAverage::sample(long iStep) 
   { 
      if (!isAtInterval(iStep)) return;

      double volume;
      volume = system().boundary().volume();
      accumulator_.sample(volume, outputFile_);

   }

   /// Output results to file after simulation is completed.
   void VolumeAverage::output() 
   {  
      // If outputFile_ was used to write block averages, close it.
      if (accumulator_.nSamplePerBlock()) {
         outputFile_.close();
      }

      fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
      ParamComposite::writeParam(outputFile_);
      outputFile_.close();

      fileMaster().openOutputFile(outputFileName(".ave"), outputFile_);
      accumulator_.output(outputFile_);
      outputFile_.close();

   }

   /*
   * Save state to binary file archive.
   */
   void VolumeAverage::save(Serializable::OArchiveType& ar)
   { ar & *this; }

   /*
   * Load state from a binary file archive.
   */
   void VolumeAverage::load(Serializable::IArchiveType& ar)
   { ar & *this; }

}

#endif 
