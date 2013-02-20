#ifndef MCMD_MD_ENERGY_OUTPUT_CPP
#define MCMD_MD_ENERGY_OUTPUT_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "MdEnergyOutput.h"
#include <mcMd/mdSimulation/md_potentials.h>
#include <util/misc/FileMaster.h>
#include <util/format/Dbl.h>
#include <util/format/Int.h>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   MdEnergyOutput::MdEnergyOutput(MdSystem& system) :
      SystemDiagnostic<MdSystem>(system)
   {  setClassName("MdEnergyOutput"); }

   /*
   * Read interval and output file name.
   */
   void MdEnergyOutput::readParameters(std::istream& in)
   {
      readInterval(in);
      readOutputFileName(in);
      fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
   }

   /*
   * Load internal state from archive.
   */
   void MdEnergyOutput::loadParameters(Serializable::IArchive &ar)
   {
      loadInterval(ar);
      loadOutputFileName(ar);
      fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
   }

   /*
   * Save internal state to archive.
   */
   void MdEnergyOutput::save(Serializable::OArchive &ar)
   { ar & *this; }

   /* 
   * Evaluate energy and print.
   */
   void MdEnergyOutput::sample(long iStep) 
   {
      if (isAtInterval(iStep)) {
	 double potential = 0.0;
         #ifndef INTER_NOPAIR
         double pair = system().pairPotential().energy();
         potential += pair;
         outputFile_ << Dbl(pair);
         #endif
         double bond = system().bondPotential().energy();
         potential += bond;
         outputFile_ << Dbl(bond);
         #ifdef INTER_ANGLE
         if (system().hasAnglePotential()) {
            double angle = system().anglePotential().energy();
            potential += angle;
            outputFile_ << Dbl(angle);
         }
         #endif
         #ifdef INTER_DIHEDRAL
         if (system().hasDihedralPotential()) {
            double dihedral = system().dihedralPotential().energy();
            potential += dihedral;
            outputFile_ << Dbl(dihedral);
         }
         #endif
         #ifdef MCMD_LINK
         if (system().hasLinkPotential()) {
            double link = system().linkPotential().energy();
            potential += link;
            outputFile_ << Dbl(link);
         }
         #endif
         #ifdef INTER_EXTERNAL
         if (system().hasExternalPotential()) {
            double external = system().externalPotential().energy();
            potential += external;
            outputFile_ << Dbl(external);
         }
         #endif
         #ifdef INTER_TETHER
         if (system().hasTetherPotential()) {
            double tether = system().tetherPotential().energy();
            potential += tether;
            outputFile_ << Dbl(tether);
         }
         #endif
         outputFile_ << Dbl(potential);
         double kinetic = system().kineticEnergy();
         outputFile_ << Dbl(kinetic);
         double total   = potential + kinetic;
         outputFile_ << Dbl(total) << std::endl;
      }
   }
 
   /* 
   * Summary
   */
   void MdEnergyOutput::output() 
   {
      // Close *.dat file
      outputFile_.close();

      // Open and write summary file
      fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
      writeParam(outputFile_);
      outputFile_ << std::endl;

      outputFile_ << "File format:" << std::endl;
      #ifndef INTER_NOPAIR
      outputFile_ << "  ";
      outputFile_ << "[pair]       ";
      #endif
      outputFile_ << "[bond]       ";
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
      if (system().hasTetherPotential()) {
         outputFile_ << "[tether]     ";
      }
      #endif
      outputFile_ << "[potential]  ";
      outputFile_ << "[kinetic]    ";
      outputFile_ << "[total]      ";
      outputFile_ << std::endl;

      outputFile_.close();
   }
}
#endif 
