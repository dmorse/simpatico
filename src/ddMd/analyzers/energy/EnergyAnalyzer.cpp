/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "EnergyAnalyzer.h"
#include <ddMd/simulation/Simulation.h>
#include <ddMd/potentials/pair/PairPotential.h>
#ifdef SIMP_BOND
#include <ddMd/potentials/bond/BondPotential.h>
#endif
#ifdef SIMP_ANGLE
#include <ddMd/potentials/angle/AnglePotential.h>
#endif
#ifdef SIMP_DIHEDRAL
#include <ddMd/potentials/dihedral/DihedralPotential.h>
#endif
#ifdef SIMP_EXTERNAL
#include <ddMd/potentials/external/ExternalPotential.h>
#endif
#include <util/accumulators/Average.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>
#include <util/mpi/MpiLoader.h>
#include <util/misc/ioUtil.h>

#include <sstream>

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   EnergyAnalyzer::EnergyAnalyzer(Simulation& simulation) 
    : Analyzer(simulation),
      totalAveragePtr_(0),
      kineticAveragePtr_(0),
      pairAveragePtr_(0),
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
   {  setClassName("EnergyAnalyzer"); }

   /*
   * Read interval and outputFileName. 
   */
   void EnergyAnalyzer::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);
      nSamplePerBlock_ = 0;
      readOptional<int>(in, "nSamplePerBlock", nSamplePerBlock_);

      Simulation& sim = simulation();
      if (sim.domain().isMaster()) {
         totalAveragePtr_ = new Average;
         totalAveragePtr_->setNSamplePerBlock(nSamplePerBlock_);
         kineticAveragePtr_ = new Average;
         kineticAveragePtr_->setNSamplePerBlock(nSamplePerBlock_);
         pairAveragePtr_ = new Average;
         pairAveragePtr_->setNSamplePerBlock(nSamplePerBlock_);
         #ifdef SIMP_BOND
         if (sim.nBondType()) {
            bondAveragePtr_ = new Average;
            bondAveragePtr_->setNSamplePerBlock(nSamplePerBlock_);
         }
         #endif
         #ifdef SIMP_ANGLE
         if (sim.nAngleType()) {
             angleAveragePtr_ = new Average;
             angleAveragePtr_->setNSamplePerBlock(nSamplePerBlock_);
         }
         #endif
         #ifdef SIMP_EXTERNAL
         if (sim.hasExternal()) {
            externalAveragePtr_ = new Average;
            externalAveragePtr_->setNSamplePerBlock(nSamplePerBlock_);
         }
         #endif
      }

      isInitialized_ = true;
   }

   /*
   * Load internal state from an archive.
   */
   void EnergyAnalyzer::loadParameters(Serializable::IArchive &ar)
   {
      loadInterval(ar);
      loadOutputFileName(ar);
      nSamplePerBlock_ = 0; // default value
      bool isRequired = false;
      loadParameter<int>(ar, "nSamplePerBlock", nSamplePerBlock_, 
                         isRequired);

      // Load Average accumulators
      Simulation& sim = simulation();
      if (sim.domain().isMaster()) {
         totalAveragePtr_ = new Average;
         ar >> *totalAveragePtr_;
         kineticAveragePtr_ = new Average;
         ar >> *kineticAveragePtr_;
         pairAveragePtr_ = new Average;
         ar >> *pairAveragePtr_;
         #ifdef SIMP_BOND
         if (sim.nBondType()) {
            bondAveragePtr_ = new Average;
            ar >> *bondAveragePtr_;
         }
         #endif
         #ifdef SIMP_ANGLE
         if (sim.nAngleType()) {
            angleAveragePtr_ = new Average;
            ar >> *angleAveragePtr_;
         }
         #endif
         #ifdef SIMP_DIHEDRAL
         if (sim.nDihedralType()) {
            dihedralAveragePtr_ = new Average;
            ar >> *dihedralAveragePtr_;
         }
         #endif
         #ifdef SIMP_EXTERNAL
         if (sim.hasExternal()) {
            externalAveragePtr_ = new Average;
            ar >> *externalAveragePtr_;
         }
         #endif
      }

      isInitialized_ = true;
   }

   /*
   * Save internal state to an archive.
   */
   void EnergyAnalyzer::save(Serializable::OArchive &ar)
   {
      saveInterval(ar);
      saveOutputFileName(ar);
      bool isActive = bool(nSamplePerBlock_);
      Parameter::saveOptional(ar, nSamplePerBlock_, isActive);

      // Save average accumulators
      Simulation& sim = simulation();
      ar << *totalAveragePtr_;
      ar << *kineticAveragePtr_;
      ar << *pairAveragePtr_;
      #ifdef SIMP_BOND
      if (sim.nBondType()) {
         ar << *bondAveragePtr_;
      }
      #endif
      #ifdef SIMP_ANGLE
      if (sim.nAngleType()) {
         ar << *angleAveragePtr_;
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (sim.nDihedralType()) {
         ar << *dihedralAveragePtr_;
      }
      #endif
      #ifdef SIMP_EXTERNAL
      if (sim.hasExternal()) {
         ar << *externalAveragePtr_;
      }
      #endif
   }
  
   /*
   * Reset nSample.
   */
   void EnergyAnalyzer::clear() 
   {  
      Simulation& sim = simulation();
      if (sim.domain().isMaster()) {
         totalAveragePtr_->clear();
         kineticAveragePtr_->clear();
         pairAveragePtr_->clear();
         #ifdef SIMP_BOND
         if (sim.nBondType()) {
            bondAveragePtr_->clear();
         }
         #endif
         #ifdef SIMP_ANGLE
         if (sim.nAngleType()) {
            angleAveragePtr_->clear();
         }
         #endif
         #ifdef SIMP_DIHEDRAL
         if (sim.nDihedralType()) {
            dihedralAveragePtr_->clear();
         }
         #endif
         #ifdef SIMP_EXTERNAL
         if (sim.hasExternal()) {
            externalAveragePtr_->clear();
         }
         #endif
      }
   }

   /*
   * Open outputfile
   */ 
   void EnergyAnalyzer::setup()
   {
      if (simulation().domain().isMaster()) {
         if (outputFile_.is_open()) {
            outputFile_.close();
         }
         std::string filename;
         filename  = outputFileName(".dat");
         simulation().fileMaster().openOutputFile(filename, outputFile_);
      }
   }

   /*
   * Output energy to file
   */
   void EnergyAnalyzer::sample(long iStep) 
   {
      if (isAtInterval(iStep))  {
         Simulation& sim = simulation();
         sim.computeKineticEnergy();
         sim.computePotentialEnergies();
         if (sim.domain().isMaster()) {
            //outputFile_ << Int(iStep, 10);
            double kinetic = sim.kineticEnergy();
            kineticAveragePtr_->sample(kinetic);
            //outputFile_ << Dbl(kinetic, 15);
            double potential = 0.0;
            double pair = sim.pairPotential().energy();
            potential += pair;
            pairAveragePtr_->sample(pair);
            //outputFile_ << Dbl(pair, 15);
            #ifdef SIMP_BOND
            if (sim.nBondType()) {
               double bond = sim.bondPotential().energy();
               potential += bond;
               bondAveragePtr_->sample(bond);
               //outputFile_ << Dbl(bond, 15);
            }
            #endif
            #ifdef SIMP_ANGLE
            if (sim.nAngleType()) {
               double angle = sim.anglePotential().energy();
               potential += angle;
               angleAveragePtr_->sample(angle);
               //outputFile_ << Dbl(angle, 15);
            }
            #endif
            #ifdef SIMP_DIHEDRAL
            if (sim.nDihedralType()) {
               double dihedral  = sim.dihedralPotential().energy();
               potential += dihedral;
               dihedralAveragePtr_->sample(dihedral);
               // outputFile_ << Dbl(dihedral, 15);
            }
            #endif
            #ifdef SIMP_EXTERNAL
            if (sim.hasExternal()) {
               double external = sim.externalPotential().energy();
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
               outputFile_ << Dbl(pairAveragePtr_->blockAverage());
               #ifdef SIMP_BOND
               if (sim.nBondType()) {
                  outputFile_ << Dbl(bondAveragePtr_->blockAverage());
               }
               #endif
               #ifdef SIMP_ANGLE
               if (sim.nAngleType()) {
                  outputFile_ << Dbl(angleAveragePtr_->blockAverage());
               }
               #endif
               #ifdef SIMP_DIHEDRAL
               if (sim.nDihedralType()) {
                  outputFile_ << Dbl(dihedralAveragePtr_->blockAverage());
               }
               #endif
               #ifdef SIMP_EXTERNAL
               if (sim.hasExternal()) {
                  outputFile_ << Dbl(externalAveragePtr_->blockAverage());
               }
               #endif
               outputFile_ << Dbl(totalAveragePtr_->blockAverage());
               outputFile_ << "\n";
            }
         }

      }
   }

   /*
   * Output results to file after simulation is completed.
   */
   void EnergyAnalyzer::output()
   {
      Simulation& sim = simulation();
      if (sim.domain().isMaster()) {

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
         ave = pairAveragePtr_->average();
         err = pairAveragePtr_->blockingError();
         outputFile_ << "Pair      " << Dbl(ave) << " +- " << Dbl(err, 9, 2) << "\n";
         #ifdef SIMP_BOND
         if (sim.nBondType()) {
            ave = bondAveragePtr_->average();
            err = bondAveragePtr_->blockingError();
            outputFile_ << "Bond      " << Dbl(ave) << " +- " << Dbl(err, 9, 2) << "\n";
         }
         #endif
         #ifdef SIMP_ANGLE
         if (sim.nAngleType()) {
            ave = angleAveragePtr_->average();
            err = angleAveragePtr_->blockingError();
            outputFile_ << "Angle     " << Dbl(ave) << " +- " << Dbl(err, 9, 2) << "\n";
         }
         #endif
         #ifdef SIMP_DIHEDRAL
         if (sim.nDihedralType()) {
            ave = dihedralAveragePtr_->average();
            err = dihedralAveragePtr_->blockingError();
            outputFile_ << "Dihedral  " << Dbl(ave) << " +- " << Dbl(err, 9, 2) << "\n";
         }
         #endif
         #ifdef SIMP_EXTERNAL
         if (sim.hasExternal()) {
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
         pairAveragePtr_->output(outputFile_);
         #ifdef SIMP_BOND
         if (sim.nBondType()) {
            outputFile_ << 
            "---------------------------------------------------------------------------------\n";
            outputFile_ << "Bond:\n\n";
            bondAveragePtr_->output(outputFile_);
         }
         #endif
         #ifdef SIMP_ANGLE
         if (sim.nAngleType()) {
            outputFile_ << 
            "---------------------------------------------------------------------------------\n";
            outputFile_ << "Angle:\n\n";
            angleAveragePtr_->output(outputFile_);
         }
         #endif
         #ifdef SIMP_DIHEDRAL
         if (sim.nDihedralType()) {
            outputFile_ << 
            "---------------------------------------------------------------------------------\n";
            outputFile_ << "Dihedral:\n\n";
            dihedralAveragePtr_->output(outputFile_);
         }
         #endif
         #ifdef SIMP_EXTERNAL
         if (sim.hasExternal()) {
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


}
