/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "MdEnergyAnalyzer.h"
#include <mcMd/mdSimulation/MdSimulation.h>
#include <mcMd/mdSimulation/MdSystem.h>

#ifndef SIMP_NOPAIR
#include <mcMd/potentials/pair/MdPairPotential.h>
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
#ifdef SIMP_COULOMB
#include <mcMd/potentials/coulomb/MdCoulombPotential.h>
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
   MdEnergyAnalyzer::MdEnergyAnalyzer(MdSystem& system) 
    : SystemAnalyzer<MdSystem>(system),
      totalAveragePtr_(0),
      kineticAveragePtr_(0),
      potentialAveragePtr_(0),
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
      #ifdef SIMP_COULOMB
      coulombRSpaceAveragePtr_(0),
      coulombKSpaceAveragePtr_(0),
      coulombAveragePtr_(0),
      coulombComponents_(false),
      #endif
      #ifdef SIMP_EXTERNAL
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
      MdSystem& sys = system();

      readInterval(in);
      readOutputFileName(in);
      nSamplePerBlock_ = 0;
      readOptional<int>(in, "nSamplePerBlock", nSamplePerBlock_);
      #ifdef SIMP_COULOMB
      if (sys.hasCoulombPotential()) {
         coulombComponents_ = false;
         readOptional<bool>(in, "coulombComponents", coulombComponents_);
      }
      #endif

      totalAveragePtr_ = new Average;
      totalAveragePtr_->setNSamplePerBlock(nSamplePerBlock_);
      kineticAveragePtr_ = new Average;
      kineticAveragePtr_->setNSamplePerBlock(nSamplePerBlock_);
      potentialAveragePtr_ = new Average;
      potentialAveragePtr_->setNSamplePerBlock(nSamplePerBlock_);
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
      #ifdef SIMP_COULOMB
      if (sys.hasCoulombPotential()) {
          if (coulombComponents_) {
             coulombRSpaceAveragePtr_ = new Average;
             coulombRSpaceAveragePtr_->setNSamplePerBlock(nSamplePerBlock_);
             coulombKSpaceAveragePtr_ = new Average;
             coulombKSpaceAveragePtr_->setNSamplePerBlock(nSamplePerBlock_);
          }
          coulombAveragePtr_ = new Average;
          coulombAveragePtr_->setNSamplePerBlock(nSamplePerBlock_);
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
   void MdEnergyAnalyzer::loadParameters(Serializable::IArchive &ar)
   {
      MdSystem& sys = system();

      // Load interval and outputFileName
      Analyzer::loadParameters(ar);

      nSamplePerBlock_ = 0; // default value
      bool isRequired = false;
      loadParameter<int>(ar, "nSamplePerBlock", nSamplePerBlock_, isRequired);
      #ifdef SIMP_COULOMB
      if (sys.hasCoulombPotential()) {
         loadParameter<bool>(ar, "coulombComponents", coulombComponents_, isRequired);
      }
      #endif

      // Load Average accumulators
      totalAveragePtr_ = new Average;
      ar >> *totalAveragePtr_;
      kineticAveragePtr_ = new Average;
      ar >> *kineticAveragePtr_;
      potentialAveragePtr_ = new Average;
      ar >> *potentialAveragePtr_;
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
      #ifdef SIMP_COULOMB
      if (sys.hasCoulombPotential()) {
         if (coulombComponents_) {
            coulombRSpaceAveragePtr_ = new Average;
            ar >> *coulombRSpaceAveragePtr_;
            coulombKSpaceAveragePtr_ = new Average;
            ar >> *coulombKSpaceAveragePtr_;
         }
         coulombAveragePtr_ = new Average;
         ar >> *coulombAveragePtr_;
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
   void MdEnergyAnalyzer::save(Serializable::OArchive &ar)
   {
      MdSystem& sys = system();

      Analyzer::save(ar);
      bool isActive = bool(nSamplePerBlock_);
      Parameter::saveOptional(ar, nSamplePerBlock_, isActive);
      #ifdef SIMP_COULOMB
      if (sys.hasCoulombPotential()) {
         isActive = coulombComponents_;
         Parameter::saveOptional(ar, coulombComponents_, isActive);
      }
      #endif

      // Save average accumulators
      ar << *totalAveragePtr_;
      ar << *kineticAveragePtr_;
      ar << *potentialAveragePtr_;
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
      #ifdef SIMP_COULOMB
      if (sys.hasCoulombPotential()) {
         if (coulombComponents_) {
            assert(coulombRSpaceAveragePtr_);
            ar << *coulombRSpaceAveragePtr_;
            assert(coulombKSpaceAveragePtr_);
            ar << *coulombKSpaceAveragePtr_;
         }
         assert(coulombAveragePtr_);
         ar << *coulombAveragePtr_;
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
   void MdEnergyAnalyzer::clear() 
   {  
      MdSystem& sys = system();

      totalAveragePtr_->clear();
      kineticAveragePtr_->clear();
      potentialAveragePtr_->clear();
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
      #ifdef SIMP_COULOMB
      if (sys.hasCoulombPotential()) {
         if (coulombComponents_) {
            assert(coulombRSpaceAveragePtr_);
            coulombRSpaceAveragePtr_->clear();
            assert(coulombKSpaceAveragePtr_);
            coulombKSpaceAveragePtr_->clear();
         }
         assert(coulombAveragePtr_);
         coulombAveragePtr_->clear();
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
   void MdEnergyAnalyzer::setup()
   {
      MdSystem& sys = system();
      Simulation& sim = sys.simulation();

      if (outputFile_.is_open()) {
         outputFile_.close();
      }

      if (nSamplePerBlock_ > 0) {

         // Open *.fmt file with format of *.dat data file
         std::string filename;
         filename  = outputFileName(".fmt");
         sim.fileMaster().openOutputFile(filename, outputFile_);

         outputFile_ << "       iStep";
         outputFile_ << "        Kinetic";
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
         #ifdef SIMP_COULOMB
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
         #ifdef SIMP_COULOMB
         if (sys.hasCoulombPotential()) {
            if (coulombComponents_) {
               // R-space contribution
               double coulombRSpace = sys.coulombPotential().rSpaceEnergy();
               assert(coulombRSpaceAveragePtr_);
               coulombRSpaceAveragePtr_->sample(coulombRSpace);
               // K-space contribution
               double coulombKSpace = sys.coulombPotential().kSpaceEnergy();
               assert(coulombKSpaceAveragePtr_);
               coulombKSpaceAveragePtr_->sample(coulombKSpace);
            }
            // Total
            double coulomb = sys.coulombPotential().energy();
            assert(coulombAveragePtr_);
            coulombAveragePtr_->sample(coulomb);
            potential += coulomb;
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

         assert(potentialAveragePtr_);
         potentialAveragePtr_->sample(potential);
         double total = kinetic + potential;
         totalAveragePtr_->sample(total);
         // outputFile_ << Dbl(total, 20)
         //             << std::endl;

         // Output block averages, if needed
         if (nSamplePerBlock_ > 0 && totalAveragePtr_->isBlockComplete()) {
            int beginStep = iStep - (nSamplePerBlock_ - 1)*interval();
            outputFile_ << Int(beginStep, 12);
            outputFile_ << Dbl(kineticAveragePtr_->blockAverage());
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
            #ifdef SIMP_COULOMB
            if (sys.hasCoulombPotential()) {
               assert(coulombAveragePtr_);
               outputFile_ << Dbl(coulombAveragePtr_->blockAverage());
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
      outputFile_ << "Kinetic    " << Dbl(ave) 
                  << " +- " << Dbl(err, 9, 2) << "\n";
      #ifndef SIMP_NOPAIR
      ave = pairAveragePtr_->average();
      err = pairAveragePtr_->blockingError();
      outputFile_ << "Pair       " << Dbl(ave) 
                  << " +- " << Dbl(err, 9, 2) << "\n";
      #endif
      #ifdef SIMP_BOND
      if (sys.hasBondPotential()) {
         ave = bondAveragePtr_->average();
         err = bondAveragePtr_->blockingError();
         outputFile_ << "Bond       " << Dbl(ave) 
                     << " +- " << Dbl(err, 9, 2) << "\n";
      }
      #endif
      #ifdef SIMP_ANGLE
      if (sys.hasAnglePotential()) {
         assert(angleAveragePtr_);
         ave = angleAveragePtr_->average();
         err = angleAveragePtr_->blockingError();
         outputFile_ << "Angle      " << Dbl(ave) 
                     << " +- " << Dbl(err, 9, 2) << "\n";
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (sys.hasDihedralPotential()) {
         assert(dihedralAveragePtr_);
         ave = dihedralAveragePtr_->average();
         err = dihedralAveragePtr_->blockingError();
         outputFile_ << "Dihedral   " << Dbl(ave) 
                     << " +- " << Dbl(err, 9, 2) << "\n";
      }
      #endif
      #ifdef SIMP_COULOMB
      if (sys.hasCoulombPotential()) {

         assert(coulombAveragePtr_);
         ave = coulombAveragePtr_->average();
         err = coulombAveragePtr_->blockingError();
         outputFile_ << "Coulomb    " << Dbl(ave) 
                     << " +- " << Dbl(err, 9, 2) << "\n";

         if (coulombComponents_) {
            assert(coulombRSpaceAveragePtr_);
            ave = coulombRSpaceAveragePtr_->average();
            err = coulombRSpaceAveragePtr_->blockingError();
            outputFile_ << "Coulomb(R) " << Dbl(ave) 
                        << " +- " << Dbl(err, 9, 2) << "\n";
            assert(coulombKSpaceAveragePtr_);
            ave = coulombKSpaceAveragePtr_->average();
            err = coulombKSpaceAveragePtr_->blockingError();
            outputFile_ << "Coulomb(K) " << Dbl(ave) << " +- " 
                        << Dbl(err, 9, 2) << "\n";
         }
      }
      #endif
      #ifdef SIMP_EXTERNAL
      if (sys.hasExternalPotential()) {
         assert(externalAveragePtr_);
         ave = externalAveragePtr_->average();
         err = externalAveragePtr_->blockingError();
         outputFile_ << "External   " << Dbl(ave) << " +- " << Dbl(err, 9, 2) << "\n";
      }
      #endif

      assert(potentialAveragePtr_);
      ave = potentialAveragePtr_->average();
      err = potentialAveragePtr_->blockingError();
      outputFile_ << "Potential  " << Dbl(ave) << " +- " << Dbl(err, 9, 2) << "\n";

      ave = totalAveragePtr_->average();
      err = totalAveragePtr_->blockingError();
      outputFile_ << "Total      " << Dbl(ave) << " +- " << Dbl(err, 9, 2) << "\n";

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
      #ifdef SIMP_COULOMB
      if (sys.hasCoulombPotential()) {
         outputFile_ << 
         "---------------------------------------------------------------------------------\n";
         outputFile_ << "Coulomb:\n\n";
         assert(coulombAveragePtr_);
         coulombAveragePtr_->output(outputFile_);
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
      outputFile_ << "Potential:\n\n";
      assert(potentialAveragePtr_);
      potentialAveragePtr_->output(outputFile_);

      outputFile_ << 
      "---------------------------------------------------------------------------------\n";
      outputFile_ << "Total:\n\n";
      assert(totalAveragePtr_);
      totalAveragePtr_->output(outputFile_);

      outputFile_.close();
   }


}
