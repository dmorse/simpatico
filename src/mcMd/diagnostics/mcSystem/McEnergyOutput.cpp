#ifndef MC_ENERGY_OUTPUT_CPP
#define MC_ENERGY_OUTPUT_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "McEnergyOutput.h"                  
#include <mcMd/mcSimulation/mc_potentials.h> // include all MC potentials
#include <mcMd/util/FileMaster.h>
#include <util/format/Dbl.h>

namespace McMd
{

   using namespace Util;

   // Constructor
   McEnergyOutput::McEnergyOutput(McSystem& system) :
      SystemDiagnostic<McSystem>(system)
   {}

   void McEnergyOutput::readParam(std::istream& in)
   {
      readInterval(in);
      readOutputFileName(in);
      fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
   }
 
   // Evaluate energy and print.
   void McEnergyOutput::sample(long iStep) 
   {
      if (isAtInterval(iStep)) {
         double energy = 0.0;
         #ifndef MCMD_NOPAIR
         double pair = system().pairPotential().energy();
         outputFile_ << Dbl(pair);
         energy += pair;
         #endif
         double bond = system().bondPotential().energy();
         outputFile_ << Dbl(bond);
         energy += bond;
         #ifdef MCMD_ANGLE
         if (system().hasAnglePotential()) {
            double angle = system().anglePotential().energy();
            outputFile_ << Dbl(angle);
            energy += angle;
         }
         #endif
         #ifdef MCMD_DIHEDRAL
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
         #ifdef MCMD_EXTERNAL
         if (system().hasExternalPotential()) {
            double external = system().externalPotential().energy();
            outputFile_ << Dbl(external);
            energy += external;
         }
         #endif
         #ifdef MCMD_TETHER
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
      #ifndef MCMD_NOPAIR
      outputFile_    << "[pair]       ";
      #endif
      outputFile_    << "[bond]       ";
      #ifdef MCMD_ANGLE
      if (system().hasAnglePotential()) {
         outputFile_ << "[angle]      ";
      }
      #endif
      #ifdef MCMD_DIHEDRAL
      if (system().hasDihedralPotential()) {
         outputFile_ << "[dihedral]   ";
      }
      #endif
      #ifdef MCMD_LINK
      if (system().hasLinkPotential()) {
         outputFile_ << "[link]       ";
      }
      #endif
      #ifdef MCMD_EXTERNAL
      if (system().hasExternalPotential()) {
         outputFile_ << "[external]   ";
      }
      #endif
      #ifdef MCMD_TETHER
      outputFile_   << "[tether]     ";
      #endif
      outputFile_    << "[potential]  ";
      outputFile_    << std::endl;

      outputFile_.close();
   }
   
}
#endif 
