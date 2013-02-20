#ifndef MCMD_MC_ENERGY_OUTPUT_CPP
#define MCMD_MC_ENERGY_OUTPUT_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "McEnergyOutput.h"                  
#include <mcMd/mcSimulation/mc_potentials.h> // include all MC potentials
#include <util/misc/FileMaster.h>
#include <util/format/Dbl.h>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   McEnergyOutput::McEnergyOutput(McSystem& system) :
      SystemDiagnostic<McSystem>(system)
   {  setClassName("McEnergyOutput"); }

   /*
   * Read file name and open output file.
   */
   void McEnergyOutput::readParameters(std::istream& in)
   {
      readInterval(in);
      readOutputFileName(in);
      fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
   }
 
   /*
   * Load state from an archive, and open output file.
   */
   void McEnergyOutput::loadParameters(Serializable::IArchive& ar)
   {  
      Diagnostic::loadParameters(ar);
      fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
   }

   /*
   * Save state to an archive.
   */
   void McEnergyOutput::save(Serializable::OArchive& ar)
   {  ar & *this; }

   /*
   * Evaluate energy and output to outputFile_.
   */
   void McEnergyOutput::sample(long iStep) 
   {
      if (isAtInterval(iStep)) {
         double energy = 0.0;
         #ifndef INTER_NOPAIR
         double pair = system().pairPotential().energy();
         outputFile_ << Dbl(pair);
         energy += pair;
         #endif
         double bond = system().bondPotential().energy();
         outputFile_ << Dbl(bond);
         energy += bond;
         #ifdef INTER_ANGLE
         if (system().hasAnglePotential()) {
            double angle = system().anglePotential().energy();
            outputFile_ << Dbl(angle);
            energy += angle;
         }
         #endif
         #ifdef INTER_DIHEDRAL
         if (system().hasDihedralPotential()) {
            double dihedral = system().dihedralPotential().energy();
            outputFile_ << Dbl(dihedral);
            energy += dihedral;
         }
         #endif
         #ifdef MCMD_LINK
         if (system().hasLinkPotential()) {
            double link = system().linkPotential().energy();
            outputFile_ << Dbl(link);
            energy += link;
         }
         #endif
         #ifdef INTER_EXTERNAL
         if (system().hasExternalPotential()) {
            double external = system().externalPotential().energy();
            outputFile_ << Dbl(external);
            energy += external;
         }
         #endif
         #ifdef INTER_TETHER
         double tether = system().tetherPotential().energy();
         outputFile_ << Dbl(tether);
         energy += tether;
         #endif
         outputFile_ << Dbl(energy) << std::endl;
      }
   }
 
   /* 
   * Summary
   */
   void McEnergyOutput::output() 
   {
      // Close *.dat file
      outputFile_.close();

      // Open and write summary file
      fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
      writeParam(outputFile_);
      outputFile_ << std::endl;
      outputFile_ << std::endl;

      outputFile_ << "File format:" << std::endl;
      #ifndef INTER_NOPAIR
      outputFile_    << "[pair]       ";
      #endif
      outputFile_    << "[bond]       ";
      #ifdef INTER_ANGLE
      if (system().hasAnglePotential()) {
         outputFile_ << "[angle]      ";
      }
      #endif
      #ifdef INTER_DIHEDRAL
      if (system().hasDihedralPotential()) {
         outputFile_ << "[dihedral]   ";
      }
      #endif
      #ifdef MCMD_LINK
      if (system().hasLinkPotential()) {
         outputFile_ << "[link]       ";
      }
      #endif
      #ifdef INTER_EXTERNAL
      if (system().hasExternalPotential()) {
         outputFile_ << "[external]   ";
      }
      #endif
      #ifdef INTER_TETHER
      outputFile_   << "[tether]     ";
      #endif
      outputFile_    << "[potential]  ";
      outputFile_    << std::endl;

      outputFile_.close();
   }
   
}
#endif 
