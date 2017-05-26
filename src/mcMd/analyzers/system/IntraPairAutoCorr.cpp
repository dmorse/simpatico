/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "IntraPairAutoCorr.h"
#include <mcMd/simulation/Simulation.h>
#include <simp/species/Species.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <util/boundary/Boundary.h>
#include <util/misc/FileMaster.h>
#include <util/archives/Serializable_includes.h>
#include <util/global.h>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   IntraPairAutoCorr::IntraPairAutoCorr(System& system) 
    : SystemAnalyzer<System>(system),
      outputFile_(),
      accumulator_(),
      data_(),
      speciesId_(-1),
      nMolecule_(-1),
      atom1Id_(-1),
      atom2Id_(-1),
      capacity_(-1),
      isInitialized_(false)
   {  setClassName("IntraPairAutoCorr"); }

   /*
   * Read parameters from file, and allocate data array.
   */
   void IntraPairAutoCorr::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);
      read<int>(in, "speciesId", speciesId_);
      read<int>(in, "atom1Id", atom1Id_);
      read<int>(in, "atom2Id", atom2Id_);
      read<int>(in, "capacity", capacity_);

      // Validate input
      if (speciesId_ < 0)       
         UTIL_THROW("Negative speciesId");
      if (atom1Id_ < 0)         
         UTIL_THROW("Negative atom1Id");
      if (atom2Id_ < 0)         
         UTIL_THROW("Negative atom2Id"); 
      if (atom2Id_ <= atom1Id_) 
         UTIL_THROW("atom2Id  <= atom1Id");
      if (capacity_ <= 0)        
         UTIL_THROW("Negative capacity");
      if (speciesId_ >= system().simulation().nSpecies()) 
         UTIL_THROW("speciesId > nSpecies");

      speciesPtr_ = &system().simulation().species(speciesId_);

      int nAtom = speciesPtr_->nAtom();
      if (atom1Id_ >= nAtom) 
         UTIL_THROW("atom1Id >= nAtom");
      if (atom2Id_ >= nAtom) 
         UTIL_THROW("atom2Id >= nAtom");

      // Maximum possible number of molecules of this species
      int speciesCapacity = speciesPtr_->capacity();

      // Allocate an array of separation Vectors
      data_.allocate(speciesCapacity); 
 
      // Allocate memory for the AutoCorrArray accumulator object
      accumulator_.setParam(speciesCapacity, capacity_);

      isInitialized_ = true;
   }

   /*
   * Load state from an archive.
   */
   void IntraPairAutoCorr::loadParameters(Serializable::IArchive& ar)
   {
      Analyzer::loadParameters(ar);
      loadParameter<int>(ar, "speciesId", speciesId_);
      loadParameter<int>(ar, "atom1Id", atom1Id_);
      loadParameter<int>(ar, "atom2Id", atom2Id_);
      loadParameter<int>(ar, "capacity", capacity_);
      ar & nMolecule_;
      ar & accumulator_;

      // Validate
      if (speciesId_ < 0) UTIL_THROW("Negative speciesId");
      if (atom1Id_ < 0) UTIL_THROW("Negative atom1Id");
      if (atom2Id_ < 0) UTIL_THROW("Negative atom2Id"); 
      if (atom2Id_ <= atom1Id_)UTIL_THROW("atom2Id  <= atom1Id");
      if (capacity_ <= 0) UTIL_THROW("Negative capacity");
      if (speciesId_ >= system().simulation().nSpecies()) {
         UTIL_THROW("speciesId > nSpecies");
      }

      speciesPtr_ = &system().simulation().species(speciesId_);
      int nAtom = speciesPtr_->nAtom();
      if (atom1Id_ >= nAtom) UTIL_THROW("atom1Id >= nAtom");
      if (atom2Id_ >= nAtom) UTIL_THROW("atom2Id >= nAtom");

      // Allocate an array of separation Vectors
      int speciesCapacity = speciesPtr_->capacity();
      data_.allocate(speciesCapacity); 
 
      isInitialized_ = true;
   }

   /*
   * Save state to archive.
   */
   void IntraPairAutoCorr::save(Serializable::OArchive& ar)
   { ar & *this; }

   /*
   * Set actual number of molecules and clear accumulator.
   */
   void IntraPairAutoCorr::setup()
   {
      if (!isInitialized_) {
         UTIL_THROW("Object not initialized.");
      }

      // Set number of molecules and clear accumulator
      nMolecule_ = system().nMolecule(speciesId_);
      accumulator_.setNEnsemble(nMolecule_);
      accumulator_.clear();
   }

   /*
   * Evaluate end-to-end vectors of all chains, add to ensemble.
   */
   void IntraPairAutoCorr::sample(long iStep) 
   { 
      if (isAtInterval(iStep))  {
         Molecule* moleculePtr;
         Vector    r1, r2, dR;
         int       i, j;

         // Validate nMolecules_
         if (nMolecule_ <= 0) {
            UTIL_THROW("nMolecule <= 0");
         }
         if (nMolecule_ != system().nMolecule(speciesId_)) {
            UTIL_THROW("Number of molecules has changed");
         }

         // Loop over molecules in species
         for (i = 0; i < nMolecule_; ++i) {
            moleculePtr = &system().molecule(speciesId_, i);

            // Build separation between atoms atom1Id and atom2Id 
            // by adding separations between consecutive atoms.
            data_[i].zero();
            for (j = atom1Id_ + 1; j <= atom2Id_; ++j) {
               r1 = moleculePtr->atom(j-1).position();
               r2 = moleculePtr->atom(j).position();
               system().boundary().distanceSq(r1, r2, dR);
               data_[i] += dR;
            }

         }
         accumulator_.sample(data_);

      } // if isAtInterval
   }

   /*
   * Output results to file after simulation is completed.
   */
   void IntraPairAutoCorr::output() 
   {  
      // Echo parameters to log file
      fileMaster().openOutputFile(outputFileName(), outputFile_);
      writeParam(outputFile_); 
      outputFile_ << std::endl;
      outputFile_ << "nMolecule       " << accumulator_.nEnsemble() 
                  << std::endl;
      outputFile_ << "bufferCapacity  " << accumulator_.bufferCapacity()
                  << std::endl;
      outputFile_ << "nSample         " << accumulator_.nSample() 
                  << std::endl;
      outputFile_ << std::endl;
      outputFile_ << "average   " << accumulator_.average() 
                  << std::endl;
      outputFile_ << std::endl;
      outputFile_ << "Format of *.dat file" << std::endl;
      outputFile_ << "[int time in samples]  [double autocorrelation function]" 
                  << std::endl;
      outputFile_ << std::endl;
      outputFile_.close();

      // Write autocorrelation function to data file
      fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      accumulator_.output(outputFile_); 
      outputFile_.close();
   }

}
