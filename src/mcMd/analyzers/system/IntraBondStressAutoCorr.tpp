#ifndef MCMD_INTRA_BOND_STRESS_AUTO_CORR_INC_H
#define MCMD_INTRA_BOND_STRESS_AUTO_CORR_INC_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "IntraBondStressAutoCorr.h"
#include <mcMd/simulation/Simulation.h>
#include <mcMd/simulation/stress.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Bond.h>
#include <mcMd/chemistry/Atom.h>
#include <simp/species/Species.h>
#include <simp/boundary/Boundary.h>
#include <util/misc/FileMaster.h>

#include <util/global.h>

namespace McMd
{

   using namespace Util;
   using namespace Simp;

   /*
   * Constructor.
   */
   template <class SystemType>
   IntraBondStressAutoCorr<SystemType>::IntraBondStressAutoCorr(SystemType& system) 
    : SystemAnalyzer<SystemType>(system),
      outputFile_(),
      accumulator_(),
      speciesId_(-1),
      nMolecule_(-1),
      nBond_(-1),
      capacity_(-1),
      isInitialized_(false)
   {}

   /*
   * Destructor.
   */
   template <class SystemType>
   IntraBondStressAutoCorr<SystemType>::~IntraBondStressAutoCorr() 
   {}

   /*
   * Read parameters from file.
   */
   template <class SystemType>
   void IntraBondStressAutoCorr<SystemType>::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);
      read(in, "speciesId", speciesId_);
      read(in, "capacity", capacity_);

      // Validate input
      if (speciesId_ < 0)       
         UTIL_THROW("Negative speciesId");
      if (speciesId_ >= system().simulation().nSpecies()) 
         UTIL_THROW("speciesId >= nSpecies");
      if (capacity_ <= 0)       
         UTIL_THROW("Negative capacity");

      speciesPtr_ = &system().simulation().species(speciesId_);
      nBond_ = speciesPtr_->nBond();
      if (nBond_ <= 0) 
         UTIL_THROW("Number of bonds per molecule <= 0"); 

      // Allocate memory
      int speciesCapacity = speciesPtr_->capacity();
      data_.allocate(speciesCapacity);
      accumulator_.setParam(speciesCapacity, capacity_);

      // Complete initialization of AutoCorrArray accumulator in
      // function setup(), when the number of molecules is known.

      isInitialized_ = true;
   }

   /*
   * Load internal state from file.
   */
   template <class SystemType> void 
   IntraBondStressAutoCorr<SystemType>::loadParameters(Serializable::IArchive &ar)
   {
      Analyzer::loadParameters(ar);
      loadParameter(ar, "speciesId", speciesId_);
      loadParameter(ar, "capacity", capacity_);
      ar & nBond_;
      ar & nMolecule_;
      ar & accumulator_;

      // Validate parameters and set speciesPtr_
      if (speciesId_ < 0) {
         UTIL_THROW("Negative speciesId");
      }
      if (speciesId_ >= system().simulation().nSpecies()) {
         UTIL_THROW("speciesId >= nSpecies");
      }
      if (capacity_ <= 0) {
         UTIL_THROW("Negative capacity");
      }
      if (nBond_ <= 0) {
         UTIL_THROW("Number of bonds per molecule <= 0");
      }
      speciesPtr_ = &system().simulation().species(speciesId_);
      if (nBond_ != speciesPtr_->nBond()) {
         UTIL_THROW("Inconsistent values for nBond");
      }

      // Allocate data_ array for internal usage.
      int speciesCapacity = speciesPtr_->capacity();
      if (speciesCapacity <= 0) {
         UTIL_THROW("Species capacity <= 0");
      }
      data_.allocate(speciesCapacity);

      isInitialized_ = true;
   }

   /*
   * Save internal state to an archive, by implicity calling serialize.
   */
   template <class SystemType> void 
   IntraBondStressAutoCorr<SystemType>::save(Serializable::OArchive &ar)
   {  ar & *this; }

   /*
   * Evaluate end-to-end vectors of all chains, add to ensemble.
   */
   template <class SystemType>
   void IntraBondStressAutoCorr<SystemType>::setup() 
   {
      // Precondition 
      UTIL_CHECK(isInitialized_);

      // Get the actual number of molecules.
      nMolecule_ = system().nMolecule(speciesId_);

      // Initialize the AutoCorrArray object
      accumulator_.setNEnsemble(nMolecule_);
      accumulator_.clear();
   }

   /*
   * Evaluate end-to-end vectors of all chains, add to ensemble.
   */
   template <class SystemType>
   void IntraBondStressAutoCorr<SystemType>::sample(long iStep) 
   { 
      UTIL_CHECK(isInitialized_);

      if (isAtInterval(iStep))  {

         Tensor    stress;
         Vector    dr, f;
         double    rsq, pressure;
         int       i, j;

         // Confirm that nMolecule has remained constant
         if (nMolecule_ != system().nMolecule(speciesId_)) {
            UTIL_THROW("Number of molecules has changed.");
         }

         // Loop over molecules
         System::ConstMoleculeIterator molIter;
         Molecule::ConstBondIterator   bondIter;
         i = 0;
         system().begin(speciesId_, molIter); 
         for ( ; molIter.notEnd(); ++ molIter) {
            data_[i].zero();

            molIter->begin(bondIter);
            for ( ; bondIter.notEnd(); ++bondIter) {
               rsq = system().boundary().distanceSq(
                                      bondIter->atom(0).position(), 
                                      bondIter->atom(1).position(), dr);
               f = dr;
               f *= system().bondPotential()
                            .forceOverR(rsq, bondIter->typeId());
               incrementPairStress(f, dr, data_[i]);
            }

            // Remove trace
            pressure  = data_[i].trace() / double(Dimension);
            for (j = 0; j < Dimension; ++j) {
               data_[i](j, j) -= pressure;
            }

            ++i;
         }

         accumulator_.sample(data_);
      } // if isAtInterval(iStep)

   }


   /*
   * Output results after simulation is completed.
   */
   template <class SystemType>
   void IntraBondStressAutoCorr<SystemType>::output() 
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
#endif 
