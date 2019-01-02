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
         Species* speciesPtr = 0;

         // Check total number of atoms
         int nFill = 0;
         for (int i = 0; i < nSpecies; ++i) {
            speciesPtr = &configuration().species(i);
            nFill += speciesPtr->nAtom()*speciesPtr->capacity();
            speciesPtr->clear();
         }
         AtomStorage* storagePtr = &configuration().atoms();
         if (storagePtr->size() != nFill) {
            std::cout << "Warning: Mismatched numbers of atoms" << std::endl;
            return false; // failure
         }

         // Set speciesId, moleculeId and atomId for each Atom
         Atom* atomPtr;
         int speciesId, moleculeId, atomId, nMolecule, nAtom;
         int id = 0;
         for (speciesId = 0; speciesId < nSpecies; ++speciesId) {
            speciesPtr = &configuration().species(speciesId);
            nMolecule  = speciesPtr->capacity();
            nAtom = speciesPtr->nAtom();
            for (moleculeId = 0; moleculeId < nMolecule; ++moleculeId) {
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

}
