#ifndef  SIMP_NOPAIR
/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "MdPairEnergyCoefficients.h"

#include <mcMd/simulation/Simulation.h>
#include <mcMd/potentials/pair/MdPairPotential.h>
#include <mcMd/neighbor/PairList.h>
#include <mcMd/neighbor/PairIterator.h>
#include <simp/species/Species.h>
#include <util/format/Dbl.h>
#include <util/format/Int.h>
#include <util/archives/Serializable_includes.h>
#include <util/misc/FileMaster.h>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   MdPairEnergyCoefficients::MdPairEnergyCoefficients(MdSystem& system)
    : SystemAnalyzer<MdSystem>(system),
    nAtomType_(system.simulation().nAtomType()),
    nSpecies_(system.simulation().nSpecies()),
    pairListPtr_(&system.pairPotential().pairList()),
    pairPotentialPtr_(&system.pairPotential()),
    boundaryPtr_(&system.boundary()),
    maxMoleculeNeighbors_(0)
   {  setClassName("MdPairEnergyCoefficients"); }

   /*
   * Destructor
   */
   MdPairEnergyCoefficients::~MdPairEnergyCoefficients()
   {}

   /*
   * Clear molecules' neighbor list arrays.
   */
   void MdPairEnergyCoefficients::clear()
   {
      int iSpecies1, iMolecule1;

      for (iSpecies1 = 0; iSpecies1 < nSpecies_; ++iSpecies1) {
         Species *species1Ptr;

         species1Ptr = &system().simulation().species(iSpecies1);
         for (iMolecule1 = 0; iMolecule1 < species1Ptr->capacity();
             ++iMolecule1) {
             moleculeNeighbors_[iSpecies1][iMolecule1].clear();
         }
      }
   }

   /*
   * Read input parameters.
   */
   void MdPairEnergyCoefficients::readParameters(std::istream& in)
   {
      readInterval(in);
      readOutputFileName(in);
      read<PairSelector>(in, "selector", selector_);
      read<int>(in, "maxMoleculeNeighbors", maxMoleculeNeighbors_);

      // Allocate
      int iSpecies;
      moleculeNeighbors_.allocate(nSpecies_);
      twoMoleculePairEnergy_.allocate(nSpecies_);
      for (iSpecies = 0; iSpecies < nSpecies_; ++iSpecies) {
         Species *speciesPtr;
         int nMolecule;
         int iMolecule;

         speciesPtr = &system().simulation().species(iSpecies); 
         nMolecule = speciesPtr->capacity();

         moleculeNeighbors_[iSpecies].allocate(nMolecule);
         twoMoleculePairEnergy_[iSpecies].allocate(nMolecule);
         for (iMolecule = 0; iMolecule < nMolecule; ++iMolecule) {
            moleculeNeighbors_[iSpecies][iMolecule].
               allocate(maxMoleculeNeighbors_);
          }
      }
      fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      isInitialized_ = true;
   }
 
   /*
   * Load state from an archive.
   */
   void MdPairEnergyCoefficients::loadParameters(Serializable::IArchive& ar)
   {
      loadInterval(ar);
      loadOutputFileName(ar);
      fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);

      loadParameter<PairSelector>(ar, "selector", selector_);
      loadParameter<int>(ar, "maxMoleculeNeighbors", maxMoleculeNeighbors_);

      int nAtomType, nSpecies;
      ar & nAtomType;
      ar & nSpecies;
     
      if (nAtomType != nAtomType_) {
         UTIL_THROW("Inconsistent values of nAtomType");
      }
      if (nSpecies != nSpecies_) {
         UTIL_THROW("Inconsistent values of nSpecies");
      }

      // Allocate moleculeNeighbors and twoMoleculePairEnergy.
      int iSpecies;
      moleculeNeighbors_.allocate(nSpecies_);
      twoMoleculePairEnergy_.allocate(nSpecies_);

      for (iSpecies = 0; iSpecies < nSpecies_; ++iSpecies) {
         Species *speciesPtr;
         int nMolecule;
         int iMolecule;

         speciesPtr = &system().simulation().species(iSpecies); 
         nMolecule = speciesPtr->capacity();

         moleculeNeighbors_[iSpecies].allocate(nMolecule);
         twoMoleculePairEnergy_[iSpecies].allocate(nMolecule);
         for (iMolecule = 0; iMolecule < nMolecule; ++iMolecule) {
            moleculeNeighbors_[iSpecies][iMolecule].
               allocate(maxMoleculeNeighbors_);
          }
      }

      ar & pairEnergyAccumulator_;
      ar & moleculePESqAccumulator_;
      ar & twoMoleculePESqAccumulator_;
      ar & pESqAccumulator_;

      isInitialized_ = true;
   }

   /*
   * Save state to archive.
   */
   void MdPairEnergyCoefficients::save(Serializable::OArchive& ar)
   { ar & *this; }

   /*
   * Evaluate contributions to accumulators.
   */
   void MdPairEnergyCoefficients::sample(long iStep) 
   {
      if (isAtInterval(iStep)) {

//       if (!system().isPairListCurrent()) {
//          system().buildPairList();
//       }

         // First clear neighbor list
         clear();

         // Loop over pairs and construct neighbor list per molecule
         PairIterator iter;
         Atom *atom0Ptr;
         Atom *atom1Ptr;

         for (pairListPtr_->begin(iter); iter.notEnd(); ++iter) {
            iter.getPair(atom0Ptr, atom1Ptr);
            if (selector_.match(*atom0Ptr, *atom1Ptr)) {
               Pair< Atom * > atomPair;
               Molecule *moleculePtr;
               Species *speciesPtr;
               int speciesId, moleculeId;

               atomPair[0] = atom0Ptr;
               atomPair[1] = atom1Ptr;

               // Store pair in neighbor list of boths atoms' molecules
               // (do double counting, since pair list counts only unique
               // pairs)
               moleculePtr = & atom0Ptr->molecule();
               speciesPtr =  & moleculePtr->species();
               speciesId = speciesPtr->id();
               moleculeId = moleculePtr->id();

               moleculeNeighbors_[speciesId][moleculeId].
                  append(atomPair);

               moleculePtr = & atom1Ptr->molecule();
               speciesPtr =  & moleculePtr->species();
               speciesId = speciesPtr->id();
               moleculeId = moleculePtr->id();

               // store atomPair in reverse order
               atomPair[0] = atom1Ptr;
               atomPair[1] = atom0Ptr;

               moleculeNeighbors_[speciesId][moleculeId].
                  append(atomPair);
            } 
         } // end pair loop

         // Now loop over all molecules
         double moleculePESq = 0.0;
         double pairEnergy=0.0;
         double twoMoleculePESq=0.0;

         double rsq;
         int iSpecies1, iMolecule1;

         for (iSpecies1 = 0; iSpecies1 < nSpecies_; ++iSpecies1) {
            Species *species1Ptr;

            species1Ptr = &system().simulation().species(iSpecies1); 
            for (iMolecule1 = 0; iMolecule1 < species1Ptr->capacity();
               ++iMolecule1) {
               int iSpecies2, iMolecule2;

               double moleculePairEnergy=0.0;
               
               // clear array for pair energy with neighboring molecules
               for (iSpecies2 = 0; iSpecies2 < nSpecies_; ++iSpecies2) {
                  Species *species2Ptr;

                  species2Ptr = &system().simulation().species(iSpecies2);
                  for (iMolecule2 = 0; iMolecule2 < species2Ptr->capacity();
                     ++iMolecule2) {
                     twoMoleculePairEnergy_[iSpecies2][iMolecule2]=0.0;
                  }
               }

               // Loop over this molecule's neighbors
               DSArray< Pair< Atom *> > *neighborsPtr;
               ConstArrayIterator< Pair< Atom *> > it;

               neighborsPtr=&moleculeNeighbors_[iSpecies1][iMolecule1];
               for (neighborsPtr->begin(it); it.notEnd(); ++it) {
                  Atom *atom0Ptr, *atom1Ptr;
                  atom0Ptr = (*it)[0];
                  atom1Ptr = (*it)[1];

                  if (selector_.match(*atom0Ptr,*atom1Ptr)) {
                     double energy;
                     Species *species2Ptr;
                     Molecule *molecule2Ptr;
                     int species2Id,molecule2Id;

                     // Load second atom's molecule and species
                     molecule2Ptr = &atom1Ptr->molecule();
                     species2Ptr = &atom1Ptr->molecule().species();
                     species2Id = species2Ptr->id();
                     molecule2Id = molecule2Ptr->id();

                     rsq = boundaryPtr_->distanceSq(atom0Ptr->position(),
                        atom1Ptr->position());
                     energy = pairPotentialPtr_->energy(rsq,
                        atom0Ptr->typeId(), atom1Ptr->typeId());

                     // increment pair energy of this molecule pair
                     twoMoleculePairEnergy_[species2Id][molecule2Id] +=
                        energy;
                  } 
               } // end neighbors loop

               // sum the molecular pair energies and their squares
               for (iSpecies2 = 0; iSpecies2 < nSpecies_; ++iSpecies2) {
                  Species *species2Ptr;

                  species2Ptr = &system().simulation().species(iSpecies2);
                  for (iMolecule2 = 0; iMolecule2 < species2Ptr->capacity();
                     ++iMolecule2) {
                     double energy = 
                        twoMoleculePairEnergy_[iSpecies2][iMolecule2];
                     moleculePairEnergy += energy;
                     twoMoleculePESq += energy*energy;
                  }

               }

               // accumulate total molecular pair energy and its square
               pairEnergy += moleculePairEnergy;
               moleculePESq += moleculePairEnergy*moleculePairEnergy;

            } // end molecule loop

         }  // end species loop

         // Sample sum of pair energies
         pairEnergyAccumulator_.sample(pairEnergy);

         // Sample sum of squares of molecular pair energy
         moleculePESqAccumulator_.sample(moleculePESq);

         // Sample sum of squares of two molecule pair energy
         twoMoleculePESqAccumulator_.sample(twoMoleculePESq);

         // Sample square of total pair energy
         pESqAccumulator_.sample(pairEnergy*pairEnergy);

         outputFile_ << Dbl(pairEnergy);
         outputFile_ << Dbl(moleculePESq);
         outputFile_ << Dbl(twoMoleculePESq);
         outputFile_ << Dbl(pairEnergy*pairEnergy);

         outputFile_ << std::endl;
      }
   }
 
   /* 
   * Summary
   */
   void MdPairEnergyCoefficients::output() 
   {
      // Close *.dat file
      outputFile_.close();

      // Open and write summary file
      fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
      writeParam(outputFile_);
      outputFile_ << std::endl;

      outputFile_ << "File format:" << std::endl;
      outputFile_ << "  ";
      outputFile_ << "[pairE]        ";
      outputFile_ << "[moPairESq]    ";
      outputFile_ << "[twoMolPairESq]";
      outputFile_ << "[pairESq]      ";
      outputFile_ << std::endl;

      outputFile_.close();

      // Write averages to separate files
      fileMaster().openOutputFile(outputFileName(".pairE.ave"),
         outputFile_);
      pairEnergyAccumulator_.output(outputFile_);
      outputFile_.close();

      fileMaster().openOutputFile(outputFileName(".molPairESq.ave"),
         outputFile_);
      moleculePESqAccumulator_.output(outputFile_);
      outputFile_.close();

      fileMaster().openOutputFile(outputFileName(".twoMolPairESq.ave"),
         outputFile_);
      twoMoleculePESqAccumulator_.output(outputFile_);
      outputFile_.close();
     
      fileMaster().openOutputFile(outputFileName(".pairESq.ave"),
         outputFile_);
      pESqAccumulator_.output(outputFile_);
      outputFile_.close();
   }

}
#endif 
