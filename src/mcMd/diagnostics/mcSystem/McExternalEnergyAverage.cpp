#ifdef INTER_EXTERNAL
#ifndef MCMD_MC_EXTERNAL_ENERGY_AVERAGE_CPP
#define MCMD_MC_EXTERNAL_ENERGY_AVERAGE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "McExternalEnergyAverage.h"        // class header

#include <mcMd/misc/FileMaster.h>  
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <mcMd/potentials/external/ExternalPotential.h>
#include <util/archives/Serializable_includes.h>


#include <cstdio> 

namespace McMd
{

   using namespace Util;

   /* 
   * Constructor.
   */
   McExternalEnergyAverage::McExternalEnergyAverage(McSystem& system)
    : SystemDiagnostic<McSystem>(system),
      outputFile_(),
      accumulator_(),
      nSamplePerBlock_(1),
      isInitialized_(false)
   {  setClassName("McExternalEnergyAverage"); }

   /*
   * Read parameters and initialize.
   */
   void McExternalEnergyAverage::readParameters(std::istream& in)
   {

      readInterval(in);
      readOutputFileName(in);
      read<int>(in,"nSamplePerBlock", nSamplePerBlock_);

      accumulator_.setNSamplePerBlock(nSamplePerBlock_);

      // If nSamplePerBlock != 0, open an output file for block averages.
      if (nSamplePerBlock_) {
         fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      }

      isInitialized_ = true;
   }

   /*
   * Clear accumulator.
   */
   void McExternalEnergyAverage::setup() 
   {  accumulator_.clear(); }
 
   /*
   * Evaluate energy per particle, and add to ensemble. 
   */
   void McExternalEnergyAverage::sample(long iStep) 
   {

      if (!isAtInterval(iStep)) return;

      double energy;
      System::MoleculeIterator molIter;
      Molecule::AtomIterator atomIter;

      energy = 0.0;

      for (int iSpec=0; iSpec < system().simulation().nSpecies(); ++iSpec){
           for (system().begin(iSpec, molIter); molIter.notEnd(); ++molIter){
                for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
                    energy += system().externalPotential().energy(atomIter->position(), atomIter->typeId());

                }
           }
      }

      accumulator_.sample(energy, outputFile_);
   }

   /*
   * Output results to file after simulation is completed.
   */
   void McExternalEnergyAverage::output() 
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
   void McExternalEnergyAverage::save(Serializable::OArchiveType& ar)
   { ar & *this; }

   /*
   * Load state from a binary file archive.
   */
   void McExternalEnergyAverage::load(Serializable::IArchiveType& ar)
   { ar & *this; }

}
#endif
#endif 
