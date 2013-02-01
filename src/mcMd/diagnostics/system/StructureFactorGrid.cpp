#ifndef MCMD_STRUCTURE_FACTOR_GRID_CPP
#define MCMD_STRUCTURE_FACTOR_GRID_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "StructureFactorGrid.h"
#include <mcMd/simulation/Simulation.h>
#include <mcMd/simulation/McMd_mpi.h>
#include <util/crystal/PointGroup.h>
#include <util/crystal/PointSymmetry.h>
#include <util/archives/Serializable_includes.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   StructureFactorGrid::StructureFactorGrid(System& system) 
    : StructureFactor(system),
      hMax_(0),
      nStar_(0),
      lattice_(Triclinic),
      isInitialized_(false),
      isFirstStep_(true)
   {  setClassName("StructureFactorGrid"); }

   /*
   * Read parameters from file, and allocate data array.
   */
   void StructureFactorGrid::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);
      read<int>(in, "nMode", nMode_);
      nAtomType_ = system().simulation().nAtomType();
      modes_.allocate(nMode_, nAtomType_);
      readDMatrix<double>(in, "modes", modes_, nMode_, nAtomType_);
      read<int>(in, "hMax", hMax_);
      read<Util::LatticeSystem>(in, "lattice", lattice_);

      // Allocate wavevectors arrays
      nWave_  = (2*hMax_ + 1)*(2*hMax_ + 1)*(2*hMax_ + 1);
      waveIntVectors_.allocate(nWave_);
      waveVectors_.allocate(nWave_);
      fourierModes_.allocate(nWave_, nMode_);
      structureFactors_.allocate(nWave_, nMode_);

      int i, j, h, k, l, m;
      IntVector g;

      // Cubic Symmetry
      if (lattice_ == Cubic) {

         nStar_ = (hMax_ + 1 )*(hMax_ + 2)*(hMax_ + 3)/6;
         starIds_.allocate(nStar_);
         starSizes_.allocate(nStar_);
   
         // Create cubic point group
         PointGroup group;
         PointSymmetry a, b, c;

         a.R(0,1) =  1;
         a.R(1,0) =  1;
         a.R(2,2) =  1;
   
         b.R(0,0) = -1;
         b.R(1,1) =  1;
         b.R(2,2) =  1;
   
         c.R(0,1) =  1;
         c.R(1,2) =  1;
         c.R(2,0) =  1;
   
         group.add(c);
         group.add(b);
         group.add(a);
         group.makeCompleteGroup();
   
         // Create grid of wavevectors
         FSArray<IntVector, 48> star;
         i = 0;
         j = 0;
         for (h = 0; h <= hMax_; ++h) {
            g[0] = h;
            for (k = 0; k <= h; ++k) {
               g[1] = k;
               for (l = 0; l <= k; ++l) {
                  g[2] = l;
                  starIds_[i] = j;
                  group.makeStar(g, star);
                  starSizes_[i] = star.size();
                  for (m = 0; m < star.size(); ++m) {
                     waveIntVectors_[j] = star[m];
                     ++j;
                  }
                  ++i;
               }
            }
         }
         if (i != nStar_) {
            UTIL_THROW("Error");
         } 
         if (j != nWave_) {
            UTIL_THROW("Error");
         } 
      } else if (lattice_ == Tetragonal) {

         nStar_ = (hMax_ + 1 )*(hMax_ + 1)*(hMax_ + 2)/2;
         starIds_.allocate(nStar_);
         starSizes_.allocate(nStar_);
         // Create tetragonal point group
         PointGroup group;
         PointSymmetry a, b, c;

         a.R(0,0) =  1;
         a.R(1,2) =  1;
         a.R(2,1) =  1;
   
         b.R(0,0) =  -1;
         b.R(1,1) =  1;
         b.R(2,2) =  1;
   
         c.R(0,0) =  1;
         c.R(1,1) =  -1;
         c.R(2,2) =  1;

         group.add(c);
         group.add(b);
         group.add(a);
         group.makeCompleteGroup();
   
         // Create grid of wavevectors
         FSArray<IntVector, 16> star;
         i = 0;
         j = 0;
         for (h = 0; h <= hMax_; ++h) {
            g[0] = h;
            for (k = 0; k <= hMax_; ++k) {
               g[1] = k;
               for (l = 0; l <= k; ++l) {
                  g[2] = l;
                  starIds_[i] = j;
                  group.makeStar(g, star);
                  starSizes_[i] = star.size();
                  for (m = 0; m < star.size(); ++m) {
                     waveIntVectors_[j] = star[m];
                     ++j;
                  }
                  ++i;
               }
            }
         }
         if (i != nStar_) {
            UTIL_THROW("Error");
         } 
         if (j != nWave_) {
            UTIL_THROW("Error");
         } 
      }

      // Clear accumulators
      for (i = 0; i < nWave_; ++i) {
         for (j = 0; j < nMode_; ++j) {
            structureFactors_(i, j) = 0.0;
         }
      }

      maximumValue_.allocate(nMode_);
      maximumWaveIntVector_.allocate(nMode_);
      maximumQ_.allocate(nMode_);
      for (int j = 0; j < nMode_; ++j) {
         maximumValue_[j].reserve(Samples);
         maximumWaveIntVector_[j].reserve(Samples);
         maximumQ_[j].reserve(Samples);
      }
     
      nSample_ = 0;

      isInitialized_ = true;
   }

   /*
   * Load state from an archive.
   */
   void StructureFactorGrid::loadParameters(Serializable::IArchive& ar)
   {
      // Load from StructureFactor::serialize
      Diagnostic::loadParameters(ar);
      ar & nAtomType_;
      loadParameter<int>(ar, "nMode", nMode_);
      loadDMatrix<double>(ar, "modes", modes_, nMode_, nAtomType_);
      ar & nWave_;
      ar & waveIntVectors_;
      ar & structureFactors_;
      ar & nSample_;

      // Load additional from StructureFactorGrid::serialize
      loadParameter<int>(ar, "hMax", hMax_);
      loadParameter<Util::LatticeSystem>(ar, "lattice", lattice_);
      ar & nStar_;
      ar & starIds_;
      ar & starSizes_;

      // Validate
      if (nWave_  != (2*hMax_ + 1)*(2*hMax_ + 1)*(2*hMax_ + 1)) {
         UTIL_THROW("Inconsistent value of nWave_");
      }
      if (nAtomType_ != system().simulation().nAtomType()) {
         UTIL_THROW("Inconsistent values of nAtomType_");
      }
      if (modes_.capacity1() != nMode_) {
         UTIL_THROW("Inconsistent capacity1 for modes array");
      }
      if (modes_.capacity2() != nAtomType_) {
         UTIL_THROW("Inconsistent capacity2 for modes array");
      }
      if (waveIntVectors_.capacity() != nWave_) {
         UTIL_THROW("Inconsistent capacity for waveIntVector");
      }

      // Allocate temporary data structures (defined in StructureFactor)
      waveVectors_.allocate(nWave_);
      fourierModes_.allocate(nWave_, nMode_);

      // Allocate maximum value history.
      maximumValue_.allocate(nMode_);
      maximumWaveIntVector_.allocate(nMode_);
      maximumQ_.allocate(nMode_);
      for (int j = 0; j < nMode_; ++j) {
         maximumValue_[j].reserve(Samples);
         maximumWaveIntVector_[j].reserve(Samples);
         maximumQ_[j].reserve(Samples);
      }

      isInitialized_ = true;
   }

   /*
   * Save state to archive.
   */
   void StructureFactorGrid::save(Serializable::OArchive& ar)
   {  ar & *this; }

   void StructureFactorGrid::setup() 
   {}

   void StructureFactorGrid::sample(long iStep)
   {
      StructureFactor::sample(iStep);

      fileMaster().openOutputFile(outputFileName(".dat"), logFile_, !isFirstStep_);
      isFirstStep_ = false;

      // Log structure factors
      double volume = system().boundary().volume();
      double norm;
      for (int i = 0; i < nStar_; ++i) {
         int size = starSizes_[i];

         int k = starIds_[i];
         for (int n = 0; n < Dimension; ++n) {
            logFile_ << Int(waveIntVectors_[k][n], 5);
         }
         logFile_ << Dbl(waveVectors_[k].abs(), 20, 8);


         for (int j = 0; j < nMode_; ++j) {
            double average = 0.0;
            double value = 0.0;
            k = starIds_[i];
            for (int m = 0; m < size; ++m) {
               norm = std::norm(fourierModes_(k, j));
               value = norm/volume;
               average += value;
               ++k;
            }
            average = average/double(size);
            logFile_ << Dbl(average, 20, 8);
         } 
         logFile_ << std::endl;
      }
      logFile_ << std::endl;
      logFile_.close();
   } 


   void StructureFactorGrid::output() 
   {
      double  value, average, size;
      int     i, j, k, m, n;

      // Echo parameters to a log file
      fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
      writeParam(outputFile_); 
      outputFile_.close();

      // Output structure factors to one file
      fileMaster().openOutputFile(outputFileName("_avg.dat"), outputFile_);

      // Loop over waves to output structure factor
      for (i = 0; i < nStar_; ++i) {
         size = starSizes_[i];

         #if 0
         // Output individual waves in star
         k = starIds_[i];
         for (m = 0; m < size; ++m) {
            for (n = 0; n < Dimension; ++n) {
               outputFile_ << Int(waveIntVectors_[k][n], 5);
            }
            outputFile_ << Dbl(waveVectors_[k].abs(), 20, 8);
            for (j = 0; j < nMode_; ++j) {
               value = structureFactors_(k, j)/double(nSample_);
               outputFile_ << Dbl(value, 20, 8);
            }
            outputFile_ << std::endl;
            ++k;
         }
         outputFile_ << std::endl;
         #endif

         k = starIds_[i];
         for (n = 0; n < Dimension; ++n) {
            outputFile_ << Int(waveIntVectors_[k][n], 5);
         }
         outputFile_ << Dbl(waveVectors_[k].abs(), 20, 8);
         for (j = 0; j < nMode_; ++j) {
            k = starIds_[i];
            average = 0.0;
            for (m = 0; m < size; ++m) {
                value = structureFactors_(k, j)/double(nSample_);
                average += value;
                ++k;
            }
            average = average/double(size);
            outputFile_ << Dbl(average, 20, 8);
         }
         outputFile_ << std::endl;

      }
      outputFile_.close();

      // Outputs history of maximum structure factors
      fileMaster().openOutputFile(outputFileName("_max.dat"), outputFile_);
      for (j = 0; j < nMode_; ++j) {
         for (int i = 0; i < nSample_; ++i) {
            outputFile_ << maximumWaveIntVector_[j][i];
            outputFile_ << Dbl(maximumQ_[j][i], 20, 8);
            outputFile_ << Dbl(maximumValue_[j][i], 20, 8);
            outputFile_ << std::endl;
         }
      }
      outputFile_.close();

   }

}
#endif
