/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McEnergyAnalyzer.h"
#include <mcMd/mcSimulation/McSimulation.h>
#include <mcMd/mcSimulation/McSystem.h>

#ifndef SIMP_NOPAIR
#include <mcMd/potentials/pair/McPairPotential.h>
#endif
#ifdef SIMP_BOND
#include <mcMd/potentials/bond/BondPotential.h>
#endif
#ifdef SIMP_ANGLE
#include <mcMd/potentials/angle/AnglePotential.h>
#endif
#ifdef SIMP_DIHEDRAL
#include <mcMd/potentials/dihedral/DihedralPotential.h>
#endif
#ifdef SIMP_EXTERNAL
#include <mcMd/potentials/external/ExternalPotential.h>
#endif

#include <util/accumulators/Average.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>
#include <util/mpi/MpiLoader.h>
#include <util/misc/ioUtil.h>

#include <sstream>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   McEnergyAnalyzer::McEnergyAnalyzer(McSystem& system) 
    : SystemAnalyzer<McSystem>(system),
      totalAveragePtr_(0),
      #ifndef SIMP_NOPAIR
      pairAveragePtr_(0),
      #endif
      #ifdef SIMP_BOND
      bondAveragePtr_(0),
      #endif
      #ifdef SIMP_ANGLE
      angleAveragePtr_(0),
      #endif
      #ifdef SIMP_DIHEDRAL
      dihedralAveragePtr_(0),
      #endif
      #ifdef SIMP_EXTERNAL
      externalAveragePtr_(0),
      #endif
      nSamplePerBlock_(0),
      isInitialized_(false)
   {  setClassName("McEnergyAnalyzer"); }

   /*
   * Read interval and outputFileName. 
   */
   void McEnergyAnalyzer::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);
      nSamplePerBlock_ = 0;
      readOptional<int>(in, "nSamplePerBlock", nSamplePerBlock_);

      McSystem& sys = system();

      totalAveragePtr_ = new Average;
      totalAveragePtr_->setNSamplePerBlock(nSamplePerBlock_);
      #ifndef SIMP_NOPAIR
      pairAveragePtr_ = new Average;
      pairAveragePtr_->setNSamplePerBlock(nSamplePerBlock_);
      #endif
      #ifdef SIMP_BOND
      if (sys.hasBondPotential()) {
         bondAveragePtr_ = new Average;
         bondAveragePtr_->setNSamplePerBlock(nSamplePerBlock_);
      }
      #endif
      #ifdef SIMP_ANGLE
      if (sys.hasAnglePotential()) {
          angleAveragePtr_ = new Average;
          angleAveragePtr_->setNSamplePerBlock(nSamplePerBlock_);
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (sys.hasDihedralPotential()) {
          dihedralAveragePtr_ = new Average;
          dihedralAveragePtr_->setNSamplePerBlock(nSamplePerBlock_);
      }
      #endif
      #ifdef SIMP_EXTERNAL
      if (sys.hasExternalPotential()) {
         externalAveragePtr_ = new Average;
         externalAveragePtr_->setNSamplePerBlock(nSamplePerBlock_);
      }
      #endif

      isInitialized_ = true;
   }

   /*
   * Load internal state from an archive.
   */
   void McEnergyAnalyzer::loadParameters(Serializable::IArchive &ar)
   {
      // Load interval and outputFileName
      Analyzer::loadParameters(ar);

      nSamplePerBlock_ = 0; // default value
      bool isRequired = false;
      loadParameter<int>(ar, "nSamplePerBlock", nSamplePerBlock_, isRequired);

      McSystem& sys = system();

      // Load Average accumulators
      totalAveragePtr_ = new Average;
      ar >> *totalAveragePtr_;
      #ifndef SIMP_NOPAIR
      pairAveragePtr_ = new Average;
      ar >> *pairAveragePtr_;
      #endif
      #ifdef SIMP_BOND
      if (sys.hasBondPotential()) {
         bondAveragePtr_ = new Average;
         ar >> *bondAveragePtr_;
      }
      #endif
      #ifdef SIMP_ANGLE
      if (sys.hasAnglePotential()) {
         angleAveragePtr_ = new Average;
         ar >> *angleAveragePtr_;
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (sys.hasDihedralPotential()) {
         dihedralAveragePtr_ = new Average;
         ar >> *dihedralAveragePtr_;
      }
      #endif
      #ifdef SIMP_EXTERNAL
      if (sys.hasExternalPotential()) {
         externalAveragePtr_ = new Average;
         ar >> *externalAveragePtr_;
      }
      #endif

      // If nSamplePerBlock != 0, open an output file for block averages.
      if (nSamplePerBlock_ != 0) {
         fileMaster().openOutputFile(Analyzer::outputFileName(".dat"), outputFile_);
      }

      isInitialized_ = true;
   }

   /*
   * Save internal state to an archive.
   */
   void McEnergyAnalyzer::save(Serializable::OArchive &ar)
   {
      Analyzer::save(ar);
      //saveInterval(ar);
      //saveOutputFileName(ar);
      bool isActive = bool(nSamplePerBlock_);
      Parameter::saveOptional(ar, nSamplePerBlock_, isActive);

      McSystem& sys = system();

      // Save average accumulators
      ar << *totalAveragePtr_;
      #ifndef SIMP_NOPAIR
      ar << *pairAveragePtr_;
      #endif
      #ifdef SIMP_BOND
      if (sys.hasBondPotential()) {
         ar << *bondAveragePtr_;
      }
      #endif
      #ifdef SIMP_ANGLE
      if (sys.hasAnglePotential()) {
         assert(angleAveragePtr_);
         ar << *angleAveragePtr_;
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (sys.hasDihedralPotential()) {
         assert(dihedralAveragePtr_);
         ar << *dihedralAveragePtr_;
      }
      #endif
      #ifdef SIMP_EXTERNAL
      if (sys.hasExternalPotential()) {
         assert(externalAveragePtr_);
         ar << *externalAveragePtr_;
      }
      #endif
   }
  
   /*
   * Reset nSample.
   */
   void McEnergyAnalyzer::clear() 
   {  
      McSystem& sys = system();

      totalAveragePtr_->clear();
      #ifndef SIMP_NOPAIR
      pairAveragePtr_->clear();
      #endif
      #ifdef SIMP_BOND
      if (sys.hasBondPotential()) {
         bondAveragePtr_->clear();
      }
      #endif
      #ifdef SIMP_ANGLE
      if (sys.hasAnglePotential()) {
         assert(angleAveragePtr_);
         angleAveragePtr_->clear();
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (sys.hasDihedralPotential()) {
         assert(dihedralAveragePtr_);
         dihedralAveragePtr_->clear();
      }
      #endif
      #ifdef SIMP_EXTERNAL
      if (sys.hasExternalPotential()) {
         assert(externalAveragePtr_);
         externalAveragePtr_->clear();
      }
      #endif
   }

   /*
   * Open outputfile
   */ 
   void McEnergyAnalyzer::setup()
   {

      if (outputFile_.is_open()) {
         outputFile_.close();
      }

      if (nSamplePerBlock_ > 0) {
         McSystem& sys = system();
         Simulation& sim = system().simulation();

         // Open *.fmt file with format of *.dat data file
         std::string filename;
         filename  = outputFileName(".fmt");
         sim.fileMaster().openOutputFile(filename, outputFile_);

         outputFile_ << "       iStep";
         #ifndef SIMP_NOPAIR
         outputFile_ << "           Pair";
         #endif
         #ifdef SIMP_BOND
         if (sys.hasBondPotential()) {
            outputFile_ << "           Bond";
         }
         #endif
         #ifdef SIMP_ANGLE
         if (sys.hasAnglePotential()) {
            outputFile_ << "          Angle";
         }
         #endif
         #ifdef SIMP_DIHEDRAL
         if (sys.hasDihedralPotential()) {
            outputFile_ << "       Dihedral";
         }
         #endif
         #if 0
         if (sys.hasCoulombPotential()) {
            outputFile_ << "        Coulomb";
         }
         #endif
         #ifdef SIMP_EXTERNAL
         if (sys.hasExternalPotential()) {
            outputFile_ << "       External";
         }
         #endif
         outputFile_ << "          Total";
         outputFile_.close();

         // Open *.dat data file, leave open for writing
         filename  = outputFileName(".dat");
         sim.fileMaster().openOutputFile(filename, outputFile_);
      }

      //std::string filename;
      //filename  = outputFileName(".dat");
      //sim.fileMaster().openOutputFile(filename, outputFile_);
   }

   /*
   * Output energy to file
   */
   void McEnergyAnalyzer::sample(long iStep) 
   {
      if (isAtInterval(iStep))  {
         McSystem& sys = system();

         //outputFile_ << Int(iStep, 10);

         double potential = 0.0;
         #ifndef SIMP_NOPAIR
         double pair = sys.pairPotential().energy();
         potential += pair;
         pairAveragePtr_->sample(pair);
         //outputFile_ << Dbl(pair, 15);
         #endif
         #ifdef SIMP_BOND
         if (sys.hasBondPotential()) {
            double bond = sys.bondPotential().energy();
            potential += bond;
            bondAveragePtr_->sample(bond);
            //outputFile_ << Dbl(bond, 15);
         }
         #endif
         #ifdef SIMP_ANGLE
         if (sys.hasAnglePotential()) {
            double angle = sys.anglePotential().energy();
            potential += angle;
            assert(angleAveragePtr_);
            angleAveragePtr_->sample(angle);
            //outputFile_ << Dbl(angle, 15);
         }
         #endif
         #ifdef SIMP_DIHEDRAL
         if (sys.hasDihedralPotential()) {
            double dihedral  = sys.dihedralPotential().energy();
            potential += dihedral;
            assert(dihedralAveragePtr_);
            dihedralAveragePtr_->sample(dihedral);
            // outputFile_ << Dbl(dihedral, 15);
         }
         #endif
         #ifdef SIMP_EXTERNAL
         if (sys.hasExternalPotential()) {
            double external = sys.externalPotential().energy();
            potential += external;
            assert(externalAveragePtr_);
            externalAveragePtr_->sample(external);
            // outputFile_ << Dbl(external, 15);
         }
         #endif
         double total = potential;
         totalAveragePtr_->sample(total);
         // outputFile_ << Dbl(total, 20)
         //             << std::endl;

         // Output block averages, if needed
         if (nSamplePerBlock_ > 0 && totalAveragePtr_->isBlockComplete()) {
            int beginStep = iStep - (nSamplePerBlock_ - 1)*interval();
            outputFile_ << Int(beginStep, 12);
            #ifndef SIMP_NOPAIR
            outputFile_ << Dbl(pairAveragePtr_->blockAverage());
            #endif
            #ifdef SIMP_BOND
            if (sys.hasBondPotential()) {
               outputFile_ << Dbl(bondAveragePtr_->blockAverage());
            }
            #endif
            #ifdef SIMP_ANGLE
            if (sys.hasAnglePotential()) {
               assert(angleAveragePtr_);
               outputFile_ << Dbl(angleAveragePtr_->blockAverage());
            }
            #endif
            #ifdef SIMP_DIHEDRAL
            if (sys.hasDihedralPotential()) {
               assert(dihedralAveragePtr_);
               outputFile_ << Dbl(dihedralAveragePtr_->blockAverage());
            }
            #endif
            #ifdef SIMP_EXTERNAL
            if (sys.hasExternalPotential()) {
               assert(externalAveragePtr_);
               outputFile_ << Dbl(externalAveragePtr_->blockAverage());
            }
            #endif
            outputFile_ << Dbl(totalAveragePtr_->blockAverage());
            outputFile_ << "\n";
         }

      }
   }

   /*
   * Output results to file after simulation is completed.
   */
   void McEnergyAnalyzer::output()
   {
      McSystem& sys = system();
      Simulation& sim = sys.simulation();

      // Close data (*.dat) file, if any
      if (outputFile_.is_open()) {
         outputFile_.close();
      }

      // Write parameter (*.prm) file
      sim.fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
      ParamComposite::writeParam(outputFile_);
      outputFile_.close();

      // Write average (*.ave) file
      sim.fileMaster().openOutputFile(outputFileName(".ave"), outputFile_);
      double ave, err;
      #ifndef SIMP_NOPAIR
      ave = pairAveragePtr_->average();
      err = pairAveragePtr_->blockingError();
      outputFile_ << "Pair      " << Dbl(ave) << " +- " << Dbl(err, 9, 2) << "\n";
      #endif
      #ifdef SIMP_BOND
      if (sys.hasBondPotential()) {
         ave = bondAveragePtr_->average();
         err = bondAveragePtr_->blockingError();
         outputFile_ << "Bond      " << Dbl(ave) << " +- " << Dbl(err, 9, 2) << "\n";
      }
      #endif
      #ifdef SIMP_ANGLE
      if (sys.hasAnglePotential()) {
         assert(angleAveragePtr_);
         ave = angleAveragePtr_->average();
         err = angleAveragePtr_->blockingError();
         outputFile_ << "Angle     " << Dbl(ave) << " +- " << Dbl(err, 9, 2) << "\n";
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (sys.hasDihedralPotential()) {
         assert(dihedralAveragePtr_);
         ave = dihedralAveragePtr_->average();
         err = dihedralAveragePtr_->blockingError();
         outputFile_ << "Dihedral  " << Dbl(ave) << " +- " << Dbl(err, 9, 2) << "\n";
      }
      #endif
      #ifdef SIMP_EXTERNAL
      if (sys.hasExternalPotential()) {
         assert(externalAveragePtr_);
         ave = externalAveragePtr_->average();
         err = externalAveragePtr_->blockingError();
         outputFile_ << "External  " << Dbl(ave) << " +- " << Dbl(err, 9, 2) << "\n";
      }
      #endif
      ave = totalAveragePtr_->average();
      err = totalAveragePtr_->blockingError();
      outputFile_ << "Total     " << Dbl(ave) << " +- " << Dbl(err, 9, 2) << "\n";
      outputFile_.close();

      // Write error analysis (*.aer) file
      sim.fileMaster().openOutputFile(outputFileName(".aer"), outputFile_);
      outputFile_ << 
      "---------------------------------------------------------------------------------\n";
      outputFile_ << "Pair:\n\n";
      #ifndef SIMP_NOPAIR
      pairAveragePtr_->output(outputFile_);
      #endif
      #ifdef SIMP_BOND
      if (sys.hasBondPotential()) {
         outputFile_ << 
         "---------------------------------------------------------------------------------\n";
         outputFile_ << "Bond:\n\n";
         bondAveragePtr_->output(outputFile_);
      }
      #endif
      #ifdef SIMP_ANGLE
      if (sys.hasAnglePotential()) {
         outputFile_ << 
         "---------------------------------------------------------------------------------\n";
         outputFile_ << "Angle:\n\n";
         assert(angleAveragePtr_);
         angleAveragePtr_->output(outputFile_);
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (sys.hasDihedralPotential()) {
         outputFile_ << 
         "---------------------------------------------------------------------------------\n";
         outputFile_ << "Dihedral:\n\n";
         assert(dihedralAveragePtr_);
         dihedralAveragePtr_->output(outputFile_);
      }
      #endif
      #ifdef SIMP_EXTERNAL
      if (sys.hasExternalPotential()) {
         outputFile_ << 
         "---------------------------------------------------------------------------------\n";
         outputFile_ << "External:\n\n";
         assert(externalAveragePtr_);
         externalAveragePtr_->output(outputFile_);
      }
      #endif
      outputFile_ << 
      "---------------------------------------------------------------------------------\n";
      outputFile_ << "Total:\n\n";
      totalAveragePtr_->output(outputFile_);
      outputFile_.close();

   }


}
