/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ConfigReader.h"
#include <mdPp/storage/Configuration.h>

#include <iostream>

namespace MdPp
{

   using namespace Util;

   /*
   * Constructor.
   */
   ConfigReader::ConfigReader()
     : configurationPtr_(0)
   {  setClassName("ConfigReader"); }

   /*
   * Constructor.
   */
   ConfigReader::ConfigReader(Configuration& configuration, 
                              bool needsAuxiliaryFile)
     : configurationPtr_(&configuration),
       needsAuxiliaryFile_(needsAuxiliaryFile)
   {  setClassName("ConfigReader"); }

   /*
   * Destructor.
   */
   ConfigReader::~ConfigReader()
   {}

   /*
   * Attempt to set atom context info, assuming ordered atom ids.
   */
   bool ConfigReader::setAtomContexts()
   {
      int nSpecies = configuration().nSpecies();
      if (nSpecies > 0) {
         SpeciesStorage* speciesPtr = 0;

         // Check total number of atoms
         int nFill = 0;
         for (int i = 0; i < nSpecies; ++i) {
            speciesPtr = &configuration().species(i);
            nFill += speciesPtr->nAtom()*speciesPtr->capacity();
         }
         AtomStorage* storagePtr = &configuration().atoms();
         if (storagePtr->size() != nFill) {
            std::cout << "Warning: Mismatched numbers of atoms" << std::endl;
            return false; // failure
         }

         // Set speciesId, moleculeId and atomId for each Atom
         Atom* atomPtr;
         int speciesId, moleculeId, atomId, nMol, nAtom;
         int id = 0;
         for (speciesId = 0; speciesId < nSpecies; ++speciesId) {
            speciesPtr = &configuration().species(speciesId);
            nMol  = speciesPtr->capacity();
            nAtom = speciesPtr->nAtom();
            for (moleculeId = 0; moleculeId < nMol; ++moleculeId) {
               for (atomId = 0; atomId < nAtom; ++atomId) {
                  atomPtr = storagePtr->ptr(id);
                  if (!atomPtr) {
                     std::cout << "Warning: missing atom" << std::endl;
                     return false;
                  }
                  atomPtr->speciesId = speciesId;
                  atomPtr->moleculeId = moleculeId;
                  atomPtr->atomId = atomId;
                  ++id;
               }
            }
         }

      }

      // Indicate successful completion
      return true;
   }

   /*
   * Add all atoms to Species containers and associated molecules.
   */
   void ConfigReader::addAtomsToSpecies()
   {
      int nSpecies = configuration().nSpecies();
      if (nSpecies > 0) {
         for (int i = 0; i < nSpecies; ++i) {
            configuration().species(i).clear();
         }
         int speciesId;
         AtomStorage::Iterator iter;
         configuration().atoms().begin(iter);
         for ( ; iter.notEnd(); ++iter) {
            speciesId = iter->speciesId;
            if (speciesId < 0) {
               UTIL_THROW("Negative speciesId");
            }
            if (speciesId >= nSpecies) {
               UTIL_THROW("speciesId >= nSpecies");
            }
            configuration().species(speciesId).addAtom(*iter);
         }
         for (int i = 0; i < nSpecies; ++i) {
            configuration().species(i).isValid();
         }
      }
   }

   #ifdef SIMP_BOND
   /*
   * Create all bonds from species templates.
   */
   void ConfigReader::makeBonds()
   {
      int nSpecies = configuration().nSpecies();
      UTIL_CHECK(nSpecies);
      GroupStorage<2>& storage = configuration().bonds();
      SpeciesStorage* speciesPtr = 0;
      int firstAtomId = 0;
      int groupId = 0;
      for (int speciesId = 0; speciesId < nSpecies; ++speciesId) {
         speciesPtr = &configuration().species(speciesId);
         int nMol = speciesPtr->capacity();
         int nAtom = speciesPtr->nAtom();
         int nGroup = speciesPtr->nBond();
         makeSpeciesGroups<2>(storage, speciesPtr->speciesBonds(),
                              nMol, nAtom, nGroup, firstAtomId, groupId);
      }
   }
   #endif

   #ifdef SIMP_ANGLE
   /*
   * Create all angles from species templates.
   */
   void ConfigReader::makeAngles()
   {
      int nSpecies = configuration().nSpecies();
      UTIL_CHECK(nSpecies);
      GroupStorage<3>& storage = configuration().angles();
      SpeciesStorage* speciesPtr = 0;
      int firstAtomId = 0;
      int groupId = 0;
      for (int speciesId = 0; speciesId < nSpecies; ++speciesId) {
         speciesPtr = &configuration().species(speciesId);
         int nMol = speciesPtr->capacity();
         int nAtom = speciesPtr->nAtom();
         int nGroup = speciesPtr->nAngle();
         makeSpeciesGroups<3>(storage, speciesPtr->speciesAngles(),
                              nMol, nAtom, nGroup, firstAtomId, groupId);
      }
   }
   #endif

   #ifdef SIMP_DIHEDRAL
   /*
   * Create all dihedrals from species templates.
   */
   void ConfigReader::makeDihedrals()
   {
      int nSpecies = configuration().nSpecies();
      UTIL_CHECK(nSpecies);
      GroupStorage<4>& storage = configuration().dihedrals();
      SpeciesStorage* speciesPtr = 0;
      int firstAtomId = 0;
      int groupId = 0;
      for (int speciesId = 0; speciesId < nSpecies; ++speciesId) {
         speciesPtr = &configuration().species(speciesId);
         int nMol = speciesPtr->capacity();
         int nAtom = speciesPtr->nAtom();
         int nGroup = speciesPtr->nDihedral();
         makeSpeciesGroups<4>(storage, speciesPtr->speciesDihedrals(),
                              nMol, nAtom, nGroup, firstAtomId, groupId);
      }
   }
   #endif

   template <int N>
   void ConfigReader::makeSpeciesGroups(
       GroupStorage<N>& storage,
       const DArray< SpeciesGroup<N> >& speciesGroups,
       int nMol, int nAtom, int nGroup, 
       int& firstAtomId, int& groupId)
   {
      Group<N>* groupPtr = 0;
      const SpeciesGroup<N>* speciesGroupPtr = 0;
      int molId, atomId, typeId, i, j;

      for (molId = 0; molId < nMol; ++molId) {
         for (i = 0; i < nGroup; ++i) {
            groupPtr = storage.newPtr();
            groupPtr->id = groupId;
            speciesGroupPtr = &speciesGroups[i]; //
            typeId = speciesGroupPtr->typeId();
            groupPtr->typeId = typeId;
            for (j = 0; j < N; ++j) {
               atomId = firstAtomId + speciesGroupPtr->atomId(j);
               groupPtr->atomIds[j] = atomId;
            }
            ++groupId;
         }
         firstAtomId += nAtom;
      }
   }

}
