/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "IntraStructureFactorGrid.h"
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
   IntraStructureFactorGrid::IntraStructureFactorGrid(System& system) 
    : IntraStructureFactor(system),
      hMax_(0),
      nStar_(0),
      lattice_(Triclinic),
      isFirstStep_(true),
      isInitialized_(false)
   {  setClassName("IntraStructureFactorGrid"); }

   /*
   * Destructor.
   */
   IntraStructureFactorGrid::~IntraStructureFactorGrid() 
   {}

   /*
   * Read parameters from file, and allocate data array.
   */
   void IntraStructureFactorGrid::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);
      read<int>(in, "speciesId", speciesId_);
      read<int>(in, "nAtomTypeIdPair", nAtomTypeIdPair_);
      nAtomType_ = system().simulation().nAtomType();
      atomTypeIdPairs_.allocate(nAtomTypeIdPair_);
      readDArray< Pair<int> >(in, "atomTypeIdPairs", atomTypeIdPairs_, 
                              nAtomTypeIdPair_);
      read<int>(in, "hMax", hMax_);
      read<Util::LatticeSystem>(in, "lattice", lattice_);

      // Allocate wavevectors arrays
      nWave_  = (2*hMax_ + 1)*(2*hMax_ + 1)*(2*hMax_ + 1);
      waveIntVectors_.allocate(nWave_);
      waveVectors_.allocate(nWave_);
      fourierModes_.allocate(nWave_, nAtomType_+1);
      structureFactors_.allocate(nWave_, nAtomTypeIdPair_);
      structureFactorDelta_.allocate(nWave_, nAtomTypeIdPair_);

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
         for (j = 0; j < nAtomTypeIdPair_; ++j) {
            structureFactors_(i, j) = 0.0;
         }
      }

      nSample_ = 0;

      isInitialized_ = true;
   }

   /*
   * Load state from an archive.
   */
   void IntraStructureFactorGrid::loadParameters(Serializable::IArchive& ar)
   {
      loadInterval(ar);
      loadOutputFileName(ar);
      loadParameter<int>(ar, "speciesId", speciesId_);
      loadParameter<int>(ar, "nAtomTypeIdPair", nAtomTypeIdPair_);
      atomTypeIdPairs_.allocate(nAtomTypeIdPair_);
      loadDArray< Pair<int> >(ar, "atomTypeIdPairs", atomTypeIdPairs_, 
                                                     nAtomTypeIdPair_);
      UTIL_CHECK(atomTypeIdPairs_.capacity() == nAtomTypeIdPair_);

      ar >> nWave_;
      //waveIntVectors_.allocate(nWave_);
      ar >> waveIntVectors_;
      //loadDArray<IntVector>(ar, "waveIntVectors", waveIntVectors_, nWave_);
      UTIL_CHECK(waveIntVectors_.capacity() == nWave_);

      ar & nAtomType_;
      UTIL_CHECK(nAtomType_ == system().simulation().nAtomType());

      // structureFactors_.allocate(nWave_, nAtomTypeIdPair_);
      ar & structureFactors_;
      UTIL_CHECK(structureFactors_.capacity1() == nWave_);
      UTIL_CHECK(structureFactors_.capacity2() == nAtomTypeIdPair_);

      //fourierModes_.allocate(nWave_, nAtomType_ + 1);
      ar & fourierModes_;
      UTIL_CHECK(fourierModes_.capacity1() == nWave_);
      UTIL_CHECK(fourierModes_.capacity2() == nAtomType_ + 1);
      ar & nSample_;

      waveVectors_.allocate(nWave_);
      structureFactorDelta_.allocate(nWave_, nAtomTypeIdPair_);

      // Load additional from IntraStructureFactorGrid::save
      loadParameter<int>(ar, "hMax", hMax_);
      UTIL_CHECK(nWave_  == (2*hMax_ + 1)*(2*hMax_ + 1)*(2*hMax_ + 1));
      loadParameter<Util::LatticeSystem>(ar, "lattice", lattice_);
      ar & nStar_;
      ar & starIds_;
      ar & starSizes_;

      isInitialized_ = true;
   }

   /*
   * Save state to archive.
   */
   void IntraStructureFactorGrid::save(Serializable::OArchive& ar)
   {  
      IntraStructureFactor::save(ar);
      ar & hMax_;
      //serializeEnum(ar, lattice_);
      ar & lattice_;
      ar & nStar_;
      ar & starIds_;
      ar & starSizes_;
   }

   void IntraStructureFactorGrid::setup() 
   {
   }

   void IntraStructureFactorGrid::sample(long iStep)
   {
      IntraStructureFactor::sample(iStep);

      // Select open mode for output files
      std::ios_base::openmode mode = std::ios_base::out;
      if (!isFirstStep_) {
        mode = std::ios_base::out | std::ios_base::app; 
      }
      fileMaster().openOutputFile(outputFileName(".log"), logFile_, mode);
      //fileMaster().openOutputFile(outputFileName(".log"), 
      //                           logFile_, !isFirstStep_);
      isFirstStep_ = false;

      // Log structure factors
      for (int i = 0; i < nStar_; ++i) {
         int size = starSizes_[i];

         int k = starIds_[i];
         for (int n = 0; n < Dimension; ++n) {
            logFile_ << Int(waveIntVectors_[k][n], 5);
         }
         logFile_ << Dbl(waveVectors_[k].abs(), 20, 8);

         for (int j = 0; j < nAtomTypeIdPair_; ++j) {
            k = starIds_[i];
            double avg = 0.0;
            for (int m = 0; m < size; ++m) {
                avg += structureFactorDelta_(k,j); 
                ++k;
            }
            avg/= size;
            logFile_ << Dbl(avg, 20, 8);
         }

         logFile_ << std::endl;
      }
      logFile_ << std::endl;
      logFile_.close();
   } 


   void IntraStructureFactorGrid::output() 
   {
      double  value, average, size;
      int     i, j, k, m, n;

      // Echo parameters to a log file
      fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
      writeParam(outputFile_); 
      outputFile_.close();

      // Output structure factors to one file
      fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);

        
      // Loop over waves to output structure factor
      for (i = 0; i < nStar_; ++i) {
         size = starSizes_[i];

         k = starIds_[i];
         for (n = 0; n < Dimension; ++n) {
            outputFile_ << Int(waveIntVectors_[k][n], 5);
         }
         outputFile_ << Dbl(waveVectors_[k].abs(), 20, 8);

         for (j = 0; j < nAtomTypeIdPair_; ++j) {
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

   }

}
