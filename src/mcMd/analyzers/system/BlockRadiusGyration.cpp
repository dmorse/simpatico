/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "BlockRadiusGyration.h"
#include <mcMd/simulation/Simulation.h>
#include <simp/species/Species.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <util/boundary/Boundary.h>
#include <util/misc/FileMaster.h>        
#include <util/archives/Serializable_includes.h>

#include <util/format/Dbl.h>

namespace McMd
{

   using namespace Util;

   /// Constructor.
   BlockRadiusGyration::BlockRadiusGyration(System& system) 
    : SystemAnalyzer<System>(system),
      outputFile_(),
      accumulators_(),
      positions_(),
      rCom_(),
      iTypeNAtom_(),
      dRSq_(),
      dRSqPair_(),
      speciesPtr_(0),
      nAtomType_(-1),
      nAtomTypePairs_(-1),
      nSamplePerBlock_(-1),
      speciesId_(-1),
      nAtom_(-1),
      isInitialized_(false)
   {  setClassName("BlockRadiusGyration"); }

   /// Read parameters from file, and allocate arrays.
   void BlockRadiusGyration::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);
      nAtomType_ = system().simulation().nAtomType();
      nAtomTypePairs_ = 0;
      int i, j;
      for (i = 0; i < nAtomType_; ++i) {
         for (j = i+1; j < nAtomType_; ++j) {
            ++nAtomTypePairs_;
         } 
      } 
      accumulators_.allocate(nAtomType_ + nAtomTypePairs_);
      
      read<int>(in,"nSamplePerBlock", nSamplePerBlock_);
      for (i = 0; i < nAtomType_+nAtomTypePairs_; ++i) {
         accumulators_[i].setNSamplePerBlock(nSamplePerBlock_);
      }
     
      // If nSamplePerBlock != 0, open an output file for block averages.
      if (accumulators_[0].nSamplePerBlock()) {
         fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      }

      read<int>(in, "speciesId", speciesId_);
      if (speciesId_ < 0) {
         UTIL_THROW("Negative speciesId");
      }

      if (speciesId_ >= system().simulation().nSpecies()) {
         UTIL_THROW("speciesId > nSpecies");
      }

      speciesPtr_ = &system().simulation().species(speciesId_);
      nAtom_ = speciesPtr_->nAtom();

      // Allocate an array of separation Vectors
      positions_.allocate(nAtom_); 
      
      // Allocate an array of center of mass position vectors
      rCom_.allocate(nAtomType_); 
      
      // Allocate an array of number of atoms in blocks of different types
      iTypeNAtom_.allocate(nAtomType_); 

      dRSq_.allocate(nAtomType_);
      dRSqPair_.allocate(nAtomTypePairs_);
 
      isInitialized_ = true;
   }

   /*
   * Load state from an archive.
   */
   void BlockRadiusGyration::loadParameters(Serializable::IArchive& ar)
   {
      // Load (everything but accumulators_)
      Analyzer::loadParameters(ar);
      loadParameter<int>(ar,"nSamplePerBlock", nSamplePerBlock_);
      loadParameter<int>(ar, "speciesId", speciesId_);
      ar & nAtom_;
      ar & nAtomType_;
      ar & nAtomTypePairs_;

      speciesPtr_ = &system().simulation().species(speciesId_);

      // Validate
      if (speciesId_ < 0) {
         UTIL_THROW("Negative speciesId");
      }
      if (speciesId_ >= system().simulation().nSpecies()) {
         UTIL_THROW("speciesId > nSpecies");
      }
      if (nAtom_ != speciesPtr_->nAtom()) {
         UTIL_THROW("Inconsistent nAtomType");
      }
      if (nAtomType_ != system().simulation().nAtomType()) {
         UTIL_THROW("Inconsistent nAtomType");
      }
      {
         int nPairTypes = 0;
         int i, j;
         for (i = 0; i < nAtomType_; ++i) {
            for (j = i+1; j < nAtomType_; ++j) {
               ++nPairTypes;
            }
         }
         if (nAtomTypePairs_ != nPairTypes) {
            UTIL_THROW("Inconsistent nAtomTypePairs");
         }
      }

      // Allocate
      positions_.allocate(nAtom_); 
      rCom_.allocate(nAtomType_); 
      iTypeNAtom_.allocate(nAtomType_); 
      dRSq_.allocate(nAtomType_);
      dRSqPair_.allocate(nAtomTypePairs_);
      accumulators_.allocate(nAtomType_ + nAtomTypePairs_);
      for (int i = 0; i < nAtomType_+nAtomTypePairs_; ++i) {
         accumulators_[i].setNSamplePerBlock(nSamplePerBlock_);
      }

      ar & accumulators_;

      // If nSamplePerBlock != 0, open an output file for block averages.
      if (accumulators_[0].nSamplePerBlock()) {
         fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      }
      isInitialized_ = true;
   }

   /*
   * Save state to archive.
   */
   void BlockRadiusGyration::save(Serializable::OArchive& ar)
   { ar & *this; }

   /*
   * Clear accumulators.
   */
   void BlockRadiusGyration::setup() 
   {  
      if (!isInitialized_) UTIL_THROW("Object is not initialized");

      for (int i = 0; i < nAtomType_ + nAtomTypePairs_; ++i) {
         accumulators_[i].clear(); 
      }
   }

   /* 
   * Evaluate end-to-end vectors of all chains, add to ensemble.
   */
   void BlockRadiusGyration::sample(long iStep) 
   { 
      if (isAtInterval(iStep))  {

         Molecule* moleculePtr;
         Vector    r1, r2, dR;
         int       i, j, k, l, m, typeId, nMolecule;

         k = 0;
         for (i = 0; i < nAtomType_; ++i) {
            dRSq_[i] = 0.0;
            iTypeNAtom_[i] = 0;
            for (j = i+1; j < nAtomType_; ++j) { 
               ++k;
               dRSqPair_[k-1] = 0.0;
            }
         }

         nMolecule = system().nMolecule(speciesId_);
         moleculePtr = &system().molecule(speciesId_,0);
         for (j = 0 ; j < nAtom_; j++) {
            typeId = moleculePtr->atom(j).typeId(); 
            iTypeNAtom_[typeId] += 1;
         }
 
         for (i = 0; i < system().nMolecule(speciesId_); i++) {
            moleculePtr = &system().molecule(speciesId_, i);

            // Initializing COM positions for different types to zero
            for (j = 0; j < nAtomType_; ++j) {
               rCom_[j].zero();
            }

            // Construct map of molecule with no periodic boundary conditions
            positions_[0] = moleculePtr->atom(0).position();
            typeId = moleculePtr->atom(0).typeId();
            rCom_[typeId] += positions_[0];
            for (j = 1 ; j < nAtom_; j++) {
               typeId = moleculePtr->atom(j).typeId(); 
               r1 = moleculePtr->atom(j-1).position();
               r2 = moleculePtr->atom(j).position();
               system().boundary().distanceSq(r1, r2, dR);
               positions_[j]  = positions_[j-1];
               positions_[j] += dR;
               rCom_[typeId] += positions_[j];
            }
            for (l = 0; l < nAtomType_; ++l) {
               rCom_[l] /= double(iTypeNAtom_[l]);
            }
            k = 0;
            for (l = 0; l < nAtomType_; ++l) {
               for (m = l+1; m < nAtomType_; ++m) {
                  ++k;
                  dR.subtract(rCom_[l], rCom_[m]);
                  dRSqPair_[k-1] += dR.square();
               }  
            }
            for (j = 0 ; j < nAtom_; j++) {
               typeId = moleculePtr->atom(j).typeId();
               dR.subtract(positions_[j], rCom_[typeId]);
               dRSq_[typeId] += dR.square();
            }
     
         }
         k = 0;
         for (i = 0; i < nAtomType_; ++i) {
            dRSq_[i] /= double(nMolecule);
            dRSq_[i] /= double(iTypeNAtom_[i]);
            accumulators_[i].sample(dRSq_[i]);
            outputFile_ << Dbl(dRSq_[i]) << "	";
            for (j = i+1; j < nAtomType_; ++j) {
               ++k;
               dRSqPair_[k-1] /= double(nMolecule);
               accumulators_[nAtomType_+k-1].sample(dRSqPair_[k-1]);
               outputFile_ << Dbl(dRSqPair_[k-1]) << "	";
            }
         }
         outputFile_ << std::endl;
      } // if isAtInterval

   }

   /*
   * Output results to file after simulation is completed.
   */
   void BlockRadiusGyration::output() 
   {
      int i, j, k; 

      // If outputFile_ was used to write block averages, close it.
      if (accumulators_[0].nSamplePerBlock()) {
         outputFile_.close();
      }

      // Write parameters to file
      fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
      ParamComposite::writeParam(outputFile_); 
      outputFile_.close();

      // Write average to file
      fileMaster().openOutputFile(outputFileName(".ave"), outputFile_);
      k = 0;
      for (i = 0; i < nAtomType_; ++i) {
         accumulators_[i].output(outputFile_);
         for (j = i+1; j < nAtomType_; ++j) {
            ++k;
            accumulators_[nAtomType_+k-1].output(outputFile_); 
            outputFile_ << std::endl;
            outputFile_ << std::endl;
         }
      }
      outputFile_.close();

   }

}
