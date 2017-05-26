/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RingRouseAutoCorr.h"
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
   RingRouseAutoCorr::RingRouseAutoCorr(System& system) 
    : SystemAnalyzer<System>(system),
      outputFile_(),
      accumulator_(),
      data_(),
      projector_(),
      speciesId_(-1),
      nMolecule_(-1),
      nAtom_(-1),
      p_(-1),
      capacity_(-1),
      isInitialized_(false)
   {  setClassName("RingRouseAutoCorr"); }

   /*
   * Destructor.
   */
   RingRouseAutoCorr::~RingRouseAutoCorr() 
   {}

   /*
   * Read parameters from file, and allocate data array.
   */
   void RingRouseAutoCorr::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);
      read<int>(in, "speciesId", speciesId_);
      read<int>(in, "p", p_);
      read<int>(in, "capacity", capacity_);

      if (speciesId_ < 0)       UTIL_THROW("Negative speciesId");
      if (p_ < 0)               UTIL_THROW("Negative mode index");
      if (capacity_ <= 0)       UTIL_THROW("Negative capacity");
      if (speciesId_ >= system().simulation().nSpecies()) 
                                UTIL_THROW("speciesId > nSpecies");

      speciesPtr_ = &system().simulation().species(speciesId_);
      int speciesCapacity = speciesPtr_->capacity();
      nAtom_ = speciesPtr_->nAtom();

      // Allocate arrays
      data_.allocate(speciesCapacity); 
      projector_.allocate(nAtom_);

      isInitialized_ = true;
   }

   /*
   * Load internal state from an archive.
   */
   void RingRouseAutoCorr::loadParameters(Serializable::IArchive &ar)
   {
      Analyzer::loadParameters(ar);
      loadParameter<int>(ar, "speciesId", speciesId_);
      loadParameter<int>(ar, "p", p_);
      loadParameter<int>(ar, "capacity", capacity_);
      ar & nAtom_;
      ar & nMolecule_;
      ar & accumulator_;
      ar & projector_;

      speciesPtr_ = &system().simulation().species(speciesId_);
      int speciesCapacity = speciesPtr_->capacity();
      data_.allocate(speciesCapacity); 
   
      // Validate input
      if (speciesId_ < 0) {
         UTIL_THROW("Negative speciesId");
      }
      if (p_ < 0) {
         UTIL_THROW("Negative mode index");
      }
      if (capacity_ <= 0) {
         UTIL_THROW("Negative capacity");
      }
      if (speciesId_ < 0) {
         UTIL_THROW("speciesId < 0");
      }
      if (speciesId_ >= system().simulation().nSpecies()) {
         UTIL_THROW("speciesId >= nSpecies");
      }
      if (nAtom_ != speciesPtr_->nAtom()) {
         UTIL_THROW("Inconsistent values of nAtom");
      }
      if (projector_.capacity() != nAtom_) {
         UTIL_THROW("Inconsistent projector capacity");
      }
      if (accumulator_.bufferCapacity() != capacity_) {
         UTIL_THROW("Inconsistent accumulator buffer capacity");
      }

      isInitialized_ = true;
   }

   /*
   * Save internal state to an archive.
   */
   void RingRouseAutoCorr::save(Serializable::OArchive &ar)
   {  ar << *this; }

   /*
   * Allocate memory, initialize accumulator, and initialize eigenvector.
   */
   void RingRouseAutoCorr::setup() 
   { 

      // Get nMolecule for this species and nAtom per molecule
      nMolecule_ = system().nMolecule(speciesId_);
      nAtom_     = speciesPtr_->nAtom();
      if (nAtom_ <= 0) UTIL_THROW("Number of atoms per molecule < 1.");

      // data_.allocate(nMolecule_); 
      // projector_.allocate(nAtom_);

      // Initialize mode projection eigenvector
      double qMode;
      if (p_ == 0) {
         for (int j = 0; j < nAtom_; ++j)
            projector_[j] = 1.0 / double(nAtom_);
      } else {
         if (p_%2 == 0) {
            qMode = 2.0*acos(-1.0)*double(p_/2)/double(nAtom_);
            for (int j = 0; j < nAtom_; ++j)
               projector_[j] = 1.0 / double(nAtom_) * cos(qMode*double(j));
         } else {
            qMode = 2.0*acos(-1.0)*double((p_+1)/2)/double(nAtom_);
            for (int j = 0; j < nAtom_; ++j)
               projector_[j] = 1.0 / double(nAtom_) * sin(qMode*double(j));
         }
      }

      // Initialize the AutoCorrArray object.
      accumulator_.setParam(nMolecule_, capacity_);

   }

   /*
   * Evaluate end-to-end vectors of all chains, add to ensemble.
   */
   void RingRouseAutoCorr::sample(long iStep) 
   { 
      if (isAtInterval(iStep))  {

         Molecule* moleculePtr;
         Vector    r1, r2, dR, atomPos;
         int       i, j;

         // Confirm that nMolecule has remained constant
         if (nMolecule_ != system().nMolecule(speciesId_)) {
            UTIL_THROW("Number of molecules has changed.");
         }

         // Loop over molecules
         for (i=0; i < nMolecule_; i++) {
            moleculePtr = &(system().molecule(speciesId_, i));

            // Retrace non-periodic shape and calculate coefficients
            atomPos = moleculePtr->atom(0).position();
            dR.multiply(atomPos, projector_[0]);
            data_[i] = dR;
            for (j = 1; j < nAtom_; j++) {
               r1 = moleculePtr->atom(j-1).position();
               r2 = moleculePtr->atom(j).position();
               system().boundary().distanceSq(r2, r1, dR);
               atomPos += dR;
               dR.multiply(atomPos, projector_[j]);
               data_[i] += dR;
            }
         }
      
         accumulator_.sample(data_);

      } // if isAtInterval

   }

   /*
   * Output results after simulation is completed.
   */
   void RingRouseAutoCorr::output() 
   {  

      // Echo parameters to analyzer log file
      fileMaster().openOutputFile(outputFileName(), outputFile_);
      writeParam(outputFile_); 
      outputFile_.close();

      // Output statistical analysis to separate data file
      fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      accumulator_.output(outputFile_); 
      outputFile_.close();

   }

}
