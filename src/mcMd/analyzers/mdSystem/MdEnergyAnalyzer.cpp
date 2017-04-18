/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "MdEnergyAnalyzer.h"
#include <mcMd/mdSimulation/MdSimulation.h>
#include <mcMd/mdSimulation/MdSystem.h>

#ifndef INTER_NOPAIR
#include <mcMd/potentials/pair/MdPairPotential.h>
#endif
#ifdef INTER_BOND
#include <mcMd/potentials/bond/BondPotential.h>
#endif
#ifdef INTER_ANGLE
#include <mcMd/potentials/angle/AnglePotential.h>
#endif
#ifdef INTER_DIHEDRAL
#include <mcMd/potentials/dihedral/DihedralPotential.h>
#endif
//#ifdef INTER_COULOMB
//#include <mcMd/potentials/coulomb/MdCoulombPotential.h>
//#endif
#ifdef INTER_EXTERNAL
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
   MdEnergyAnalyzer::MdEnergyAnalyzer(MdSystem& system) 
    : SystemAnalyzer<MdSystem>(system),
      totalAveragePtr_(0),
      kineticAveragePtr_(0),
      #ifndef INTER_NOPAIR
      pairAveragePtr_(0),
      #endif
      #ifdef INTER_BOND
      bondAveragePtr_(0),
      #endif
      #ifdef INTER_ANGLE
      angleAveragePtr_(0),
      #endif
      #ifdef INTER_DIHEDRAL
      dihedralAveragePtr_(0),
      #endif
      #ifdef INTER_EXTERNAL
      externalAveragePtr_(0),
      #endif
      nSamplePerBlock_(0),
      isInitialized_(false)
   {  setClassName("MdEnergyAnalyzer"); }

   /*
   * Read interval and outputFileName. 
   */
   void MdEnergyAnalyzer::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);
      nSamplePerBlock_ = 0;
      readOptional<int>(in, "nSamplePerBlock", nSamplePerBlock_);

      MdSystem& sys = system();

      totalAveragePtr_ = new Average;
      totalAveragePtr_->setNSamplePerBlock(nSamplePerBlock_);
      kineticAveragePtr_ = new Average;
      kineticAveragePtr_->setNSamplePerBlock(nSamplePerBlock_);
      #ifndef INTER_NOPAIR
      pairAveragePtr_ = new Average;
      pairAveragePtr_->setNSamplePerBlock(nSamplePerBlock_);
      #endif
      #ifdef INTER_BOND
      if (sys.hasBondPotential()) {
         bondAveragePtr_ = new Average;
         bondAveragePtr_->setNSamplePerBlock(nSamplePerBlock_);
      }
      #endif
      #ifdef INTER_ANGLE
      if (sys.hasAnglePotential()) {
          angleAveragePtr_ = new Average;
          angleAveragePtr_->setNSamplePerBlock(nSamplePerBlock_);
      }
      #endif
      #ifdef INTER_DIHEDRAL
      if (sys.hasDihedralPotential()) {
          dihedralAveragePtr_ = new Average;
          dihedralAveragePtr_->setNSamplePerBlock(nSamplePerBlock_);
      }
      #endif
      #ifdef INTER_EXTERNAL
      if (sys.hasExternal()) {
         externalAveragePtr_ = new Average;
         externalAveragePtr_->setNSamplePerBlock(nSamplePerBlock_);
      }
      #endif

      isInitialized_ = true;
   }

   /*
   * Load internal state from an archive.
   */
   void MdEnergyAnalyzer::loadParameters(Serializable::IArchive &ar)
   {
      // Load interval and outputFileName
      Analyzer::loadParameters(ar);

      nSamplePerBlock_ = 0; // default value
      bool isRequired = false;
      loadParameter<int>(ar, "nSamplePerBlock", nSamplePerBlock_, isRequired);

      MdSystem& sys = system();

      // Load Average accumulators
      totalAveragePtr_ = new Average;
      ar >> *totalAveragePtr_;
      kineticAveragePtr_ = new Average;
      ar >> *kineticAveragePtr_;
      UTIL_CHECK(kineticAveragePtr_->nSamplePerBlock() == nSamplePerBlock_);
      #ifndef INTER_NOPAIR
      pairAveragePtr_ = new Average;
      ar >> *pairAveragePtr_;
      #endif
      #ifdef INTER_BOND
      if (sys.hasBondPotential()) {
         bondAveragePtr_ = new Average;
         ar >> *bondAveragePtr_;
      }
      #endif
      #ifdef INTER_ANGLE
      if (sys.hasAnglePotential()) {
         angleAveragePtr_ = new Average;
         ar >> *angleAveragePtr_;
      }
      #endif
      #ifdef INTER_DIHEDRAL
      if (sys.hasDihedralPotential()) {
         dihedralAveragePtr_ = new Average;
         ar >> *dihedralAveragePtr_;
      }
      #endif
      #ifdef INTER_EXTERNAL
      if (sys.hasExternal()) {
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
   void MdEnergyAnalyzer::save(Serializable::OArchive &ar)
   {
      Analyzer::save(ar);
      bool isActive = bool(nSamplePerBlock_);
      Parameter::saveOptional(ar, nSamplePerBlock_, isActive);

      MdSystem& sys = system();

      // Save average accumulators
      ar << *totalAveragePtr_;
      ar << *kineticAveragePtr_;
      #ifndef INTER_NOPAIR
      ar << *pairAveragePtr_;
      #endif
      #ifdef INTER_BOND
      if (sys.hasBondPotential()) {
         ar << *bondAveragePtr_;
      }
      #endif
      #ifdef INTER_ANGLE
      if (sys.hasAnglePotential()) {
         ar << *angleAveragePtr_;
      }
      #endif
      #ifdef INTER_DIHEDRAL
      if (sys.hasDihedralPotential()) {
         ar << *dihedralAveragePtr_;
      }
      #endif
      #ifdef INTER_EXTERNAL
      if (sys.hasExternal()) {
         ar << *externalAveragePtr_;
      }
      #endif
   }
  
   /*
   * Reset nSample.
   */
   void MdEnergyAnalyzer::clear() 
   {  
      MdSystem& sys = system();

      totalAveragePtr_->clear();
      kineticAveragePtr_->clear();
      #ifndef INTER_NOPAIR
      pairAveragePtr_->clear();
      #endif
      #ifdef INTER_BOND
      if (sys.hasBondPotential()) {
         bondAveragePtr_->clear();
      }
      #endif
      #ifdef INTER_ANGLE
      if (sys.hasAnglePotential()) {
         angleAveragePtr_->clear();
      }
      #endif
      #ifdef INTER_DIHEDRAL
      if (sys.hasDihedralPotential()) {
         dihedralAveragePtr_->clear();
      }
      #endif
      #ifdef INTER_EXTERNAL
      if (sys.hasExternal()) {
         externalAveragePtr_->clear();
      }
      #endif
   }

   /*
   * Open outputfile
   */ 
   void MdEnergyAnalyzer::setup()
   {
      Simulation& sim = system().simulation();

      if (outputFile_.is_open()) {
         outputFile_.close();
      }
      std::string filename;
      filename  = outputFileName(".dat");
      sim.fileMaster().openOutputFile(filename, outputFile_);
   }

   /*
   * Output energy to file
   */
   void MdEnergyAnalyzer::sample(long iStep) 
   {
      if (isAtInterval(iStep))  {
         MdSystem& sys = system();

         //outputFile_ << Int(iStep, 10);
         double kinetic = sys.kineticEnergy();
         kineticAveragePtr_->sample(kinetic);
         //outputFile_ << Dbl(kinetic, 15);

         double potential = 0.0;
         #ifndef INTER_NOPAIR
         double pair = sys.pairPotential().energy();
         potential += pair;
         pairAveragePtr_->sample(pair);
         //outputFile_ << Dbl(pair, 15);
         #endif
         #ifdef INTER_BOND
         if (sys.hasBondPotential()) {
            double bond = sys.bondPotential().energy();
            potential += bond;
            bondAveragePtr_->sample(bond);
            //outputFile_ << Dbl(bond, 15);
         }
         #endif
         #ifdef INTER_ANGLE
         if (sys.hasAnglePotential()) {
            double angle = sys.anglePotential().energy();
            potential += angle;
            angleAveragePtr_->sample(angle);
            //outputFile_ << Dbl(angle, 15);
         }
         #endif
         #ifdef INTER_DIHEDRAL
         if (sys.hasDihedralPotential()) {
            double dihedral  = sys.dihedralPotential().energy();
            potential += dihedral;
            dihedralAveragePtr_->sample(dihedral);
            // outputFile_ << Dbl(dihedral, 15);
         }
         #endif
         #ifdef INTER_EXTERNAL
         if (sys.hasExternal()) {
            double external = sys.externalPotential().energy();
            potential += external;
            externalAveragePtr_->sample(external);
            // outputFile_ << Dbl(external, 15);
         }
         #endif
         double total = kinetic + potential;
         totalAveragePtr_->sample(total);
         // outputFile_ << Dbl(total, 20)
         //             << std::endl;

         // Output block averages, if needed
         if (nSamplePerBlock_ > 0 && totalAveragePtr_->isBlockComplete()) {
            int beginStep = iStep - (nSamplePerBlock_ - 1)*interval();
            outputFile_ << Int(beginStep, 12);
            outputFile_ << Dbl(kineticAveragePtr_->blockAverage());
            #ifndef INTER_NOPAIR
            outputFile_ << Dbl(pairAveragePtr_->blockAverage());
            #endif
            #ifdef INTER_BOND
            if (sys.hasBondPotential()) {
               outputFile_ << Dbl(bondAveragePtr_->blockAverage());
            }
            #endif
            #ifdef INTER_ANGLE
            if (sys.hasAnglePotential()) {
               outputFile_ << Dbl(angleAveragePtr_->blockAverage());
            }
            #endif
            #ifdef INTER_DIHEDRAL
            if (sys.hasDihedralPotential()) {
               outputFile_ << Dbl(dihedralAveragePtr_->blockAverage());
            }
            #endif
            #ifdef INTER_EXTERNAL
            if (sys.hasExternal()) {
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
   void MdEnergyAnalyzer::output()
   {
      MdSystem& sys = system();
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
      ave = kineticAveragePtr_->average();
      err = kineticAveragePtr_->blockingError();
      outputFile_ << "Kinetic   " << Dbl(ave) << " +- " << Dbl(err, 9, 2) << "\n";
      #ifndef INTER_NOPAIR
      ave = pairAveragePtr_->average();
      err = pairAveragePtr_->blockingError();
      outputFile_ << "Pair      " << Dbl(ave) << " +- " << Dbl(err, 9, 2) << "\n";
      #endif
      #ifdef INTER_BOND
      if (sys.hasBondPotential()) {
         ave = bondAveragePtr_->average();
         err = bondAveragePtr_->blockingError();
         outputFile_ << "Bond      " << Dbl(ave) << " +- " << Dbl(err, 9, 2) << "\n";
      }
      #endif
      #ifdef INTER_ANGLE
      if (sys.hasAnglePotential()) {
         ave = angleAveragePtr_->average();
         err = angleAveragePtr_->blockingError();
         outputFile_ << "Angle     " << Dbl(ave) << " +- " << Dbl(err, 9, 2) << "\n";
      }
      #endif
      #ifdef INTER_DIHEDRAL
      if (sys.hasDihedralPotential()) {
         ave = dihedralAveragePtr_->average();
         err = dihedralAveragePtr_->blockingError();
         outputFile_ << "Dihedral  " << Dbl(ave) << " +- " << Dbl(err, 9, 2) << "\n";
      }
      #endif
      #ifdef INTER_EXTERNAL
      if (sys.hasExternal()) {
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
      outputFile_ << "Kinetic:\n\n";
      kineticAveragePtr_->output(outputFile_);
      outputFile_ << 
      "---------------------------------------------------------------------------------\n";
      outputFile_ << "Pair:\n\n";
      #ifndef INTER_NOPAIR
      pairAveragePtr_->output(outputFile_);
      #endif
      #ifdef INTER_BOND
      if (sys.hasBondPotential()) {
         outputFile_ << 
         "---------------------------------------------------------------------------------\n";
         outputFile_ << "Bond:\n\n";
         bondAveragePtr_->output(outputFile_);
      }
      #endif
      #ifdef INTER_ANGLE
      if (sys.hasAnglePotential()) {
         outputFile_ << 
         "---------------------------------------------------------------------------------\n";
         outputFile_ << "Angle:\n\n";
         angleAveragePtr_->output(outputFile_);
      }
      #endif
      #ifdef INTER_DIHEDRAL
      if (sys.hasDihedralPotential()) {
         outputFile_ << 
         "---------------------------------------------------------------------------------\n";
         outputFile_ << "Dihedral:\n\n";
         dihedralAveragePtr_->output(outputFile_);
      }
      #endif
      #ifdef INTER_EXTERNAL
      if (sys.hasExternal()) {
         outputFile_ << 
         "---------------------------------------------------------------------------------\n";
         outputFile_ << "External:\n\n";
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
