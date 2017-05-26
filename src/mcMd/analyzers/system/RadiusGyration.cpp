/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RadiusGyration.h"
#include <mcMd/simulation/Simulation.h>
#include <simp/species/Species.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <util/boundary/Boundary.h>
#include <util/misc/FileMaster.h>        
#include <util/archives/Serializable_includes.h>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   RadiusGyration::RadiusGyration(System& system) 
    : SystemAnalyzer<System>(system),
      isInitialized_(false)
   {  setClassName("RadiusGyration"); }

   /*
   * Read parameters, allocate memory, and initialize accumulator.
   */
   void RadiusGyration::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);
      read<int>(in,"nSamplePerBlock", nSamplePerBlock_);
      read<int>(in, "speciesId", speciesId_);

      if (speciesId_ < 0) {
         UTIL_THROW("Negative speciesId");
      }
      if (speciesId_ >= system().simulation().nSpecies()) {
         UTIL_THROW("speciesId > nSpecies");
      }
      speciesPtr_ = &system().simulation().species(speciesId_);
      nAtom_ = speciesPtr_->nAtom();

      positions_.allocate(nAtom_); 
      accumulator_.setNSamplePerBlock(nSamplePerBlock_);
      accumulator_.clear();

      // Open output file for block averages, if nSamplePerBlock != 0.
      if (accumulator_.nSamplePerBlock()) {
         fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      }
      isInitialized_ = true;
   }

   /*
   * Load state from an archive.
   */
   void RadiusGyration::loadParameters(Serializable::IArchive& ar)
   {
      Analyzer::loadParameters(ar);
      loadParameter<int>(ar,"nSamplePerBlock", nSamplePerBlock_);
      loadParameter<int>(ar, "speciesId", speciesId_);
      ar & nAtom_;

      // Validate
      if (speciesId_ < 0) {
         UTIL_THROW("Negative speciesId");
      }
      if (speciesId_ >= system().simulation().nSpecies()) {
         UTIL_THROW("speciesId > nSpecies");
      }
      speciesPtr_ = &system().simulation().species(speciesId_);
      if (nAtom_ != speciesPtr_->nAtom()) {
         UTIL_THROW("Inconsistent values for nAtom");
      }

      ar & accumulator_;
      ar & positions_;

      // Open output file for block averages, if nSamplePerBlock != 0.
      if (nSamplePerBlock_) {
         fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      }
      isInitialized_ = true;
   }

   /*
   * Save state to archive.
   */
   void RadiusGyration::save(Serializable::OArchive& ar)
   { ar & *this; }

   /*
   * Clear accumulator (public method).
   */
   void RadiusGyration::setup() 
   {
      if (!isInitialized_) {
         UTIL_THROW("Object is not initialized");
      }  
      accumulator_.clear(); 
   }

   /* 
   * Evaluate squared radii of gyration of all chains, add to ensemble.
   */
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

            // Construct unwrapped map of molecule (no periodic b.c.'s)
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

            // Calculate dRSq
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
   * Output final results to file, after simulation is completed.
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

}
