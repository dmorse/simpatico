#ifndef MCMD_MC_MU_EXCHANGE_CPP
#define MCMD_MC_MU_EXCHANGE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McMuExchange.h"                            // class header
#include <util/misc/FileMaster.h>
#include <util/archives/Serializable_includes.h>

#include <mcMd/mcSimulation/McSystem.h>
#include <mcMd/potentials/pair/McPairPotential.h>
#include <mcMd/neighbor/CellList.h>

#include <mcMd/species/Species.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <util/boundary/Boundary.h>
#include <util/ensembles/EnergyEnsemble.h>
#include <util/space/Vector.h>
#include <util/global.h>

#include <math.h>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   McMuExchange::McMuExchange(McSystem& system)
    : SystemAnalyzer<McSystem>(system),
      simulationPtr_(&system.simulation()),
      boundaryPtr_(&system.boundary()),
      outputFile_(),
      accumulators_(),
      newTypeIds_(),
      flipAtomIds_(),
      speciesId_(-1),
      nAtom_(-1),
      isInitialized_(false)
   {}


   /*
   * Read parameters and initialize.
   */
   void McMuExchange::readParameters(std::istream& in)
   {

      readInterval(in);
      readOutputFileName(in);
      read<int>(in, "speciesId", speciesId_);
      if (speciesId_ < 0) {
         UTIL_THROW("Negative speciesId");
      }
      if (speciesId_ >= system().simulation().nSpecies()) {
         UTIL_THROW("speciesId >= nSpecies");
      }
      Species* speciesPtr;
      speciesPtr = &(system().simulation().species(speciesId_));
      nAtom_ = speciesPtr->nAtom();
      newTypeIds_.allocate(nAtom_);
      readDArray<int>(in, "newTypeIds", newTypeIds_, nAtom_);

      flipAtomIds_.allocate(nAtom_);
      for (int i = 0; i < nAtom_; ++i) {
         if (newTypeIds_[i] != speciesPtr->atomTypeId(i)) {
            flipAtomIds_.append(i);
         }
      }
      accumulators_.allocate(speciesPtr->capacity());

      isInitialized_ = true;
   }

   /*
   * Clear accumulators.
   */
   void McMuExchange::setup()
   {
      nMolecule_ = system().nMolecule(speciesId_);
      for (int iMol = 0; iMol < nMolecule_; ++iMol) {
         accumulators_[iMol].clear(); 
      }
   }

   /*
   * Evaluate change in energy, add Boltzmann factor to accumulator.
   */
   void McMuExchange::sample(long iStep)
   {
      if (isAtInterval(iStep))  {

         // Precondition
         if (!system().energyEnsemble().isIsothermal()) {
            UTIL_THROW("EnergyEnsemble is not isothermal");
         }

         McPairPotential& potential = system().pairPotential();
         const CellList& cellList = potential.cellList();
         System::MoleculeIterator molIter;
         Atom* ptr0 = 0;    // Pointer to first atom in molecule
         Atom* ptr1 = 0;    // Pointer to flipped atom
         Atom* ptr2 = 0;    // Pointer to neighbor
         Mask* maskPtr = 0; // Mask of flipped atom
         double beta = system().energyEnsemble().beta();
         double rsq, dE, boltzmann;
         int j, k, nNeighbor;
         int i1, i2, id1, id2, t1, t2, t1New, iMol;

         // Loop over molecules in species
         iMol = 0;
         system().begin(speciesId_, molIter); 
         for ( ; molIter.notEnd(); ++molIter) {
            dE = 0.0;
            ptr0 = &molIter->atom(0);

            // Loop over flipped atoms
            for (j=0; flipAtomIds_.size(); ++j) {
               i1 = flipAtomIds_[j];
               t1New = newTypeIds_[i1];
               ptr1 = &molIter->atom(i1);
               id1 = ptr1->id();
               t1 = ptr1->typeId();
               maskPtr = &(ptr1->mask());

               // Loop over neighboring atoms
               cellList.getNeighbors(ptr1->position(), neighbors_);
               nNeighbor = neighbors_.size();
               for (k = 0; k < nNeighbor; ++k) {
                  ptr2 = neighbors_[k];
                  id2 = ptr2->id();

                  // Check if atoms are identical
                  if (id2 != id1) {
          
                     // Check if pair is masked
                     if (!maskPtr->isMasked(*ptr2)) {

                        rsq = boundary().distanceSq(ptr1->position(), 
                                                    ptr2->position());
                        t2 = ptr2->typeId();
                        if (&(ptr1->molecule()) != &(ptr2->molecule())) {

                           // Intermolecular atom pair 
                           dE -= potential.energy(rsq, t1, t2);
                           dE += potential.energy(rsq, t1New, t2);

                        } else {

                           // Intramolecular atom pair 
                           if (id2 > id1) {
                              dE -= potential.energy(rsq, t1, t2);
                              i2 = (int)(ptr2 - ptr0);
                              t2 = newTypeIds_[i2];
                              dE += potential.energy(rsq, t1New, t2);
                           }

                        }
                     }
                  }
               } // end loop over neighbors
            } // end loop over flipped atoms
            boltzmann = exp(-beta*dE);
            accumulators_[iMol].sample(boltzmann);
            ++iMol;
         } // end loop over molecules
      }
   }

   /*
   * Output results to file after simulation is completed.
   */
   void McMuExchange::output()
   {

      fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
      writeParam(outputFile_);
      outputFile_.close();

      fileMaster().openOutputFile(outputFileName(".ave"), outputFile_);
      for (int iMol=0; iMol < nMolecule_; ++iMol) {
         accumulators_[iMol].output(outputFile_);
      }
      outputFile_.close();
   }

   /*
   * Save state to binary file archive.
   */
   void McMuExchange::save(Serializable::OArchive& ar)
   {  ar & *this; }

   /*
   * Load state from a binary file archive.
   */
   void McMuExchange::load(Serializable::IArchive& ar)
   { ar & *this; }

}
#endif
