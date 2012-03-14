#ifndef INTRA_BOND_TENSOR_AUTO_CORR_INC_H
#define INTRA_BOND_TENSOR_AUTO_CORR_INC_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "IntraBondTensorAutoCorr.h"
#include <mcMd/simulation/Simulation.h>
#include <mcMd/species/Species.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Bond.h>
#include <mcMd/chemistry/Atom.h>
#include <util/boundary/Boundary.h>
#include <mcMd/util/FileMaster.h>

#include <util/global.h>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   template <class SystemType>
   IntraBondTensorAutoCorr<SystemType>::IntraBondTensorAutoCorr(SystemType& system) 
    : SystemDiagnostic<SystemType>(system),
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
   IntraBondTensorAutoCorr<SystemType>::~IntraBondTensorAutoCorr() 
   {}

   /*
   * Read parameters from file.
   */
   template <class SystemType>
   void IntraBondTensorAutoCorr<SystemType>::readParam(std::istream& in) 
   {

      // Read interval and parameters for AutoCorrArray
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
      if (nBond_ <= 0) UTIL_THROW("Number of bonds per molecule <= 0");

      // Allocate memory
      int speciesCapacity = speciesPtr_->capacity();
      if (speciesCapacity <= 0) UTIL_THROW("Species capacity <= 0");
      data_.allocate(speciesCapacity);
      accumulator_.setParam(speciesCapacity, capacity_);

      // Delay initialization of AutoCorrArray until first call
      // of sample(), when the number of molecules is known.

      isInitialized_ = true;
   }

   /*
   * Evaluate end-to-end vectors of all chains, add to ensemble.
   */
   template <class SystemType>
   void IntraBondTensorAutoCorr<SystemType>::initialize() 
   { 

      if (!isInitialized_) {
         UTIL_THROW("Object is not intitialized");
      }

      // Get number of molecules and number of atoms per molecule.
      nMolecule_ = system().nMolecule(speciesId_);
      if (nMolecule_ <= 0) UTIL_THROW("nMolecule <= 0");

      // Initialize the AutoCorrArray object
      accumulator_.setNEnsemble(nMolecule_);
   }

   /*
   * Evaluate end-to-end vectors of all chains, add to ensemble.
   */
   template <class SystemType>
   void IntraBondTensorAutoCorr<SystemType>::sample(long iStep) 
   { 
      if (isAtInterval(iStep))  {

         Tensor  t;
         Vector  dr, u;
         double  trace;
         int     i, j;

         // Confirm that nMolecule has remained constant
         if (nMolecule_ != system().nMolecule(speciesId_)) {
            UTIL_THROW("Number of molecules has changed.");
         }

         Boundary& boundary = system().boundary();
         System::ConstMoleculeIterator molIter;
         Molecule::ConstBondIterator   bondIter;

         // Loop over molecules
         i = 0;
         system().begin(speciesId_, molIter); 
         for ( ; molIter.notEnd(); ++ molIter) {
            data_[i].zero();

            // Loop over bonds
            molIter->begin(bondIter);
            for ( ; bondIter.notEnd(); ++bondIter) {
               boundary.distanceSq(bondIter->atom(0).position(), 
                                   bondIter->atom(1).position(), dr);
               u.versor(dr);
               t.dyad(u, u);
               data_[i] += t;
            }

            // Remove trace
            trace  = data_[i].trace()/double(Dimension);
            for (j=0; j < Dimension; ++j) {
               data_[i](j, j) -= trace;
            }

            ++i;
         }

         accumulator_.sample(data_);

      } // Is at interval
   }

   /// Output results after simulation is completed.
   template <class SystemType>
   void IntraBondTensorAutoCorr<SystemType>::output() 
   {  

      // Echo parameters to diagnostic log file
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
