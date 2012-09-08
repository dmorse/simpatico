#ifndef MCMD_RADIUS_GYRATION_CPP
#define MCMD_RADIUS_GYRATION_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "RadiusGyration.h"
#include <mcMd/simulation/Simulation.h>
#include <mcMd/species/Species.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <util/boundary/Boundary.h>
#include <mcMd/util/FileMaster.h>        
#include <util/archives/Serializable_includes.h>

namespace McMd
{

   using namespace Util;

   /// Constructor.
   RadiusGyration::RadiusGyration(System& system) 
    : SystemDiagnostic<System>(system),
      isInitialized_(false)
   {  setClassName("RadiusGyration"); }

   /// Read parameters from file, and allocate data array.
   void RadiusGyration::readParameters(std::istream& in) 
   {

      readInterval(in);
      readOutputFileName(in);
      read<int>(in,"nSamplePerBlock", nSamplePerBlock_);

      accumulator_.setNSamplePerBlock(nSamplePerBlock_);
      accumulator_.clear();

      // If nSamplePerBlock != 0, open an output file for block averages.
      if (accumulator_.nSamplePerBlock()) {
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

      isInitialized_ = true;
   }

   /*
   * Clear accumulator.
   */
   void RadiusGyration::setup() 
   {
      if (!isInitialized_) {
         UTIL_THROW("Object is not initialized");
      }  
      accumulator_.clear(); 
   }
 
   /// Evaluate end-to-end vectors of all chains, add to ensemble.
   void RadiusGyration::sample(long iStep) 
   { 
      if (isAtInterval(iStep))  {

         Molecule* moleculePtr;
         Vector    r1, r2, dR, Rcm;
         double    dRSq;
         int       i, j, nMolecule;

         dRSq = 0.0;
         nMolecule = system().nMolecule(speciesId_);
         for (i = 0; i < system().nMolecule(speciesId_); i++) {
            moleculePtr = &system().molecule(speciesId_, i);

            // Construct map of molecule with no periodic boundary conditions
            positions_[0] = moleculePtr->atom(0).position();
            Rcm = positions_[0];
            for (j = 1 ; j < nAtom_; j++) {
               r1 = moleculePtr->atom(j-1).position();
               r2 = moleculePtr->atom(j).position();
               system().boundary().distanceSq(r1, r2, dR);
               positions_[j]  = positions_[j-1];
               positions_[j] += dR;
               Rcm += positions_[j];
            }
            Rcm /= double(nAtom_);

            for (j = 0 ; j < nAtom_; j++) {
               dR.subtract(positions_[j], Rcm);
               dRSq += dR.square();
            }
     
         }
         dRSq /= double(nMolecule);
         dRSq /= double(nAtom_);
      
         accumulator_.sample(dRSq, outputFile_);

      } 

   }

   /*
   * Output results to file after simulation is completed.
   */
   void RadiusGyration::output() 
   { 
      // If outputFile_ was used to write block averages, close it.
      if (accumulator_.nSamplePerBlock()) {
         outputFile_.close();
      }

      // Write parameters to file
      fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
      ParamComposite::writeParam(outputFile_); 
      outputFile_.close();

      // Write average to file
      fileMaster().openOutputFile(outputFileName(".ave"), outputFile_);
      accumulator_.output(outputFile_); 
      outputFile_.close();
   }

   /*
   * Save state to binary file archive.
   */
   void RadiusGyration::save(Serializable::OArchiveType& ar)
   { ar & *this; }

   /*
   * Load state from a binary file archive.
   */
   void RadiusGyration::load(Serializable::IArchiveType& ar)
   { ar & *this; }
}
#endif 
