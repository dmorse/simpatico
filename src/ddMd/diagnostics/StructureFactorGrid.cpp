#ifndef DDMD_STRUCTURE_FACTOR_GRID_CPP
#define DDMD_STRUCTURE_FACTOR_GRID_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "StructureFactorGrid.h"
#include <ddMd/simulation/Simulation.h>
#include <util/crystal/PointGroup.h>
#include <util/crystal/PointSymmetry.h>
#include <util/archives/Serializable_includes.h>
#include <util/mpi/MpiLoader.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

#include <iostream>
#include <fstream>
namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   StructureFactorGrid::StructureFactorGrid(Simulation& simulation) 
    : StructureFactor(simulation),
      hMax_(0),
      nStar_(0),
      lattice_(Triclinic)
   {  setClassName("StructureFactorGrid"); }

   /*
   * Read parameters from file, and allocate data array.
   */
   void StructureFactorGrid::readParameters(std::istream& in) 
   {
      nAtomType_ = simulation().nAtomType();

      // Read parameters
      readInterval(in);
      readOutputFileName(in);
      read<int>(in, "nMode", nMode_);
      modes_.allocate(nMode_, nAtomType_);
      readDMatrix<double>(in, "modes", modes_, nMode_, nAtomType_);
      read<int>(in, "hMax", hMax_);
      read<LatticeSystem>(in, "lattice", lattice_);

      // Allocate wavevectors arrays
      nWave_     = (2*hMax_ +1 )*(2*hMax_ + 1)*(2*hMax_ + 1);
      waveIntVectors_.allocate(nWave_);
      waveVectors_.allocate(nWave_);
      fourierModes_.allocate(nWave_, nMode_);
      totalFourierModes_.allocate(nWave_, nMode_);

      int i, j, h, k, l, m;
      IntVector g;

      // Cubic Symmetry
      if (lattice_ == Cubic) {
         nStar_ = (hMax_ +1 )*(hMax_ + 2)*(hMax_ + 3)/6;
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
      } else if (lattice_ == Orthorhombic) {

         nStar_ = (hMax_ + 1 )*(hMax_ + 1)*(hMax_ + 1);
         starIds_.allocate(nStar_);
         starSizes_.allocate(nStar_);
         // Create tetragonal point group
         PointGroup group;
         PointSymmetry a, b, c;

         a.R(0,0) =  -1;
         a.R(1,1) =  1;
         a.R(2,2) =  1;

         b.R(0,0) =  1;
         b.R(1,1) =  -1;
         b.R(2,2) =  1;

         c.R(0,0) =  1;
         c.R(1,1) =  1;
         c.R(2,2) =  -1;

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

      if (simulation().domain().isMaster()) {

         structureFactors_.allocate(nWave_, nMode_);
         int i, j;
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

      }
      nSample_ = 0;

      isInitialized_ = true;
   }

   /*
   * Load internal state from an archive.
   */
   void StructureFactorGrid::loadParameters(Serializable::IArchive &ar)
   {
      nAtomType_ = simulation().nAtomType();

      // Load parameter file parameters
      loadInterval(ar);
      loadOutputFileName(ar);
      loadParameter<int>(ar, "nMode", nMode_);
      modes_.allocate(nMode_, nAtomType_);
      loadDMatrix<double>(ar, "modes", modes_, nMode_, nAtomType_);
      loadParameter<int>(ar, "hMax", hMax_);
      loadParameter<LatticeSystem>(ar, "lattice", lattice_);

      // Load and broadcast other distributed members
      MpiLoader<Serializable::IArchive> loader(*this, ar);
      loader.load(nWave_);
      waveIntVectors_.allocate(nWave_);
      loader.load(waveIntVectors_);
      loader.load(nStar_);
      starIds_.allocate(nStar_);
      loader.load(starIds_);
      starSizes_.allocate(nStar_);
      loader.load(starSizes_);
      loader.load(nSample_);

      if (simulation().domain().isMaster()) {

         structureFactors_.allocate(nWave_, nMode_);
         ar >> structureFactors_;

         maximumValue_.allocate(nMode_);
         maximumWaveIntVector_.allocate(nMode_);
         maximumQ_.allocate(nMode_);
         for (int j = 0; j < nMode_; ++j) {
            ar >> maximumValue_[j];
            ar >> maximumWaveIntVector_[j];
            ar >> maximumQ_[j];
         }
      }

      // Allocate work space
      waveVectors_.allocate(nWave_);
      fourierModes_.allocate(nWave_, nMode_);
      totalFourierModes_.allocate(nWave_, nMode_);

      isInitialized_ = true;
   }

   /*
   * Save internal state to an archive.
   */
   void StructureFactorGrid::save(Serializable::OArchive &ar)
   {
      saveInterval(ar);
      saveOutputFileName(ar);
      ar << nMode_;
      ar << modes_;
      ar << hMax_;
      ar << lattice_;

      ar << nWave_;
      ar << waveIntVectors_;
      ar << nStar_;
      ar << starIds_;
      ar << starSizes_;
      ar << nSample_;

      ar << structureFactors_;
      for (int j = 0; j < nMode_; ++j) {
         ar << maximumValue_[j];
         ar << maximumWaveIntVector_[j];
         ar << maximumQ_[j];
      }
   }

   void StructureFactorGrid::output()
   {
      if (simulation().domain().isMaster()) {
            
         double  value, average, size;
         int     i, j, k, m, n;

         // Write parameters to a *.prm file
         simulation().fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
         writeParam(outputFile_);
         outputFile_.close();

         // Output all structure factors to one file
         simulation().fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
         for (i = 0; i < nStar_; ++i) {
            size = starSizes_[i];
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

         // Output history of maximum structure factors to *_max.dat file
         simulation().fileMaster().openOutputFile(outputFileName("_max.dat"), 
                                                  outputFile_);
         for (j = 0; j < nMode_; ++j) {
            for (i = 0; i < nSample_; ++i) {
               for (n = 0; n < Dimension; ++n) {
                  outputFile_ << Int(maximumWaveIntVector_[j][i][n], 5);
               }
               outputFile_ << Dbl(maximumQ_[j][i], 20, 8);
               outputFile_ << Dbl(maximumValue_[j][i], 20, 8);
               outputFile_ << std::endl;
            }
         }
         outputFile_.close();
      
      }
   }

}
#endif
