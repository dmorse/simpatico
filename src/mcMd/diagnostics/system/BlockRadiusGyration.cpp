#ifndef MCMD_BLOCK_RADIUS_GYRATION_CPP
#define MCMD_BLOCK_RADIUS_GYRATION_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "BlockRadiusGyration.h"
#include <mcMd/simulation/Simulation.h>
#include <mcMd/species/Species.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <util/boundary/Boundary.h>
#include <mcMd/util/FileMaster.h>        
#include <util/archives/Serializable_includes.h>

#include <util/format/Dbl.h>

namespace McMd
{

   using namespace Util;

   /// Constructor.
   BlockRadiusGyration::BlockRadiusGyration(System& system) 
    : SystemDiagnostic<System>(system),
      outputFile_(),
      accumulators_(),
      positions_(),
      rCom_(),
      iTypeNAtom_(),
      speciesPtr_(0),
      nAtomType_(-1),
      nAtomTypePairs_(-1),
      nSamplePerBlock_(-1),
      speciesId_(-1),
      nAtom_(-1),
      isInitialized_(false)
   {  setClassName("BlockRadiusGyration"); }

   /// Read parameters from file, and allocate data array.
   void BlockRadiusGyration::readParam(std::istream& in) 
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

      isInitialized_ = true;
   }

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
         DArray<double> dRSq;      //to store dRSq of different atomTypes
         DArray<double> dRSqPair;  //to store dRSq of different atomTypePairs
         int       i, j, k, l, m, typeId, nMolecule;

         dRSq.allocate(nAtomType_);
         dRSqPair.allocate(nAtomTypePairs_);
         k = 0;
         for (i = 0; i < nAtomType_; ++i) {
            dRSq[i] = 0.0;
            iTypeNAtom_[i] = 0;
            for (j = i+1; j < nAtomType_; ++j) { 
               ++k;
               dRSqPair[k-1] = 0.0;
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

            //Initializing the center of mass positions for different block types to be zero
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
                  dRSqPair[k-1] += dR.square();
               }  
            }
            for (j = 0 ; j < nAtom_; j++) {
               typeId = moleculePtr->atom(j).typeId();
               dR.subtract(positions_[j], rCom_[typeId]);
               dRSq[typeId] += dR.square();
            }
     
         }
         k = 0;
         for (i = 0; i < nAtomType_; ++i) {
            dRSq[i] /= double(nMolecule);
            dRSq[i] /= double(iTypeNAtom_[i]);
            accumulators_[i].sample(dRSq[i]);
            outputFile_ << Dbl(dRSq[i]) << "	";
            for (j = i+1; j < nAtomType_; ++j) {
               ++k;
               dRSqPair[k-1] /= double(nMolecule);
               accumulators_[nAtomType_+k-1].sample(dRSqPair[k-1]);
               outputFile_ << Dbl(dRSqPair[k-1]) << "	";
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

   /*
   * Save state to binary file archive.
   */
   void BlockRadiusGyration::save(Serializable::OArchiveType& ar)
   { ar & *this; }

   /*
   * Load state from a binary file archive.
   */
   void BlockRadiusGyration::load(Serializable::IArchiveType& ar)
   { ar & *this; }

}
#endif 
