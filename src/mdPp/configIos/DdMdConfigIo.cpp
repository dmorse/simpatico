#ifndef MDPP_DDMD_CONFIG_IO_CPP
#define MDPP_DDMD_CONFIG_IO_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "DdMdConfigIo.h"

#include <mdPp/chemistry/Atom.h>
#include <mdPp/chemistry/Group.h>
#include <mdPp/chemistry/Species.h>
//#include <mdPp/chemistry/MaskPolicy.h>
#include <mdPp/storage/Storage.h>

#include <util/space/Vector.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

namespace MdPp
{

   using namespace Util;

   /*
   * Constructor.
   */
   DdMdConfigIo::DdMdConfigIo(Storage& storage, bool hasMolecules)
    : ConfigIo(storage),
      hasMolecules_(hasMolecules)
   {  setClassName("DdMdConfigIo"); }

   /*
   * Read a configuration file.
   */
   void DdMdConfigIo::readConfig(std::ifstream& file)
   {
      // Precondition
      if (!file.is_open()) {  
            UTIL_THROW("Error: File is not open"); 
      }

      // Read and broadcast boundary
      file >> Label("BOUNDARY");
      file >> storage().boundary();

      // Read and distribute atoms

      // Read atoms
      Atom* atomPtr;
      int atomCapacity = storage().atoms().capacity(); // Maximum allowed id + 1
      int nAtom;          
      file >> Label("ATOMS");
      file >> Label("nAtom") >> nAtom;
      for (int i = 0; i < nAtom; ++i) {

         // Get pointer to new atom 
         atomPtr = storage().atoms().newPtr();
 
         file >> atomPtr->id >> atomPtr->typeId;
         if (atomPtr->id < 0 || atomPtr->id >= atomCapacity) {
            UTIL_THROW("Invalid atom id");
         }
         if (hasMolecules_) {
            file >> atomPtr->speciesId 
                 >> atomPtr->moleculeId 
                 >> atomPtr->atomId;
            if (atomPtr->speciesId < 0) {
               UTIL_THROW("Invalid species id");
            }
            if (atomPtr->moleculeId < 0) {
               UTIL_THROW("Invalid molecule id");
            }
            if (atomPtr->atomId < 0) {
               UTIL_THROW("Invalid atom id");
            }
         }
         file >> atomPtr->position;
         file >> atomPtr->velocity;

         storage().atoms().add();
      }

      // Read Covalent Groups
      #ifdef INTER_BOND
      if (storage().bonds().capacity()) {
         readGroups(file, "BONDS", "nBond", storage().bonds());
         //if (maskPolicy == MaskBonded) {
         //   setAtomMasks();
         //}
      }
      #endif

      // Optionally add atoms to species
      if (storage().nSpecies() > 0) {
         if (!hasMolecules_) {
            UTIL_THROW("No atom context info in chosen ConfigIo");
         }
         int speciesId;
         AtomStorage::Iterator iter;
         storage().atoms().begin(iter);
         for ( ; iter.notEnd(); ++iter) {
            speciesId = iter->speciesId;
            storage().species(speciesId).addAtom(*iter);
         }
         for (int i = 0; storage().nSpecies(); ++i) {
            storage().species(i).isValid();
         }
      }
   }

   /* 
   * Write the configuration file.
   */
   void DdMdConfigIo::writeConfig(std::ofstream& file)
   {
      // Precondition
      if (!file.is_open()) {  
         UTIL_THROW("Error: File is not open"); 
      }

      // Write Boundary dimensions
      file << "BOUNDARY" << std::endl << std::endl;
      file << storage().boundary() << std::endl;
      file << std::endl;

      // Atoms
      int nAtom = storage().atoms().size();
      file << "ATOMS" << std::endl;
      file << "nAtom" << Int(nAtom, 10) << std::endl;

      // Write atoms
      Vector r;
      AtomStorage::Iterator iter;
      storage().atoms().begin(iter);
      for (; iter.notEnd(); ++iter) {
         file << Int(iter->id, 10) 
              << Int(iter->typeId, 6);
         r = iter->position;
         if (hasMolecules_) {
            file << Int(iter->speciesId, 6) 
                 << Int(iter->moleculeId, 10)
                 << Int(iter->atomId, 6);
         }
         file << "\n" << r 
              << "\n" << iter->velocity << "\n";
      }

      // Write the groups
      #ifdef INTER_BOND
      if (storage().bonds().capacity()) {
         writeGroups(file, "BONDS", "nBond", storage().bonds());
      }
      #endif

   }
 
}
#endif
