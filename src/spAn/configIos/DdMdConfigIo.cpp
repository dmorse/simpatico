#ifndef SPAN_DDMD_CONFIG_IO_CPP
#define SPAN_DDMD_CONFIG_IO_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "DdMdConfigIo.h"

#include <spAn/chemistry/Atom.h>
#include <spAn/chemistry/Group.h>
#include <spAn/chemistry/Species.h>
//#include <spAn/chemistry/MaskPolicy.h>
#include <spAn/storage/Configuration.h>

#include <util/space/Vector.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

namespace SpAn
{

   using namespace Util;

   /*
   * Constructor.
   */
   DdMdConfigIo::DdMdConfigIo(Configuration& configuration, bool hasMolecules)
    : ConfigIo(configuration),
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
      file >> configuration().boundary();

      // Read and distribute atoms

      // Read atoms
      Atom* atomPtr;
      int atomCapacity = configuration().atoms().capacity(); // Maximum allowed id + 1
      int nAtom;          
      file >> Label("ATOMS");
      file >> Label("nAtom") >> nAtom;
      for (int i = 0; i < nAtom; ++i) {

         // Get pointer to new atom 
         atomPtr = configuration().atoms().newPtr();
 
         file >> atomPtr->id;
         if (atomPtr->id < 0) {
            std::cout << "atom id =" << atomPtr->id << std::endl;
            UTIL_THROW("Negative atom id");
         }
         if (atomPtr->id >= atomCapacity) {
            std::cout << "atom id      =" << atomPtr->id << std::endl;
            std::cout << "atomCapacity =" << atomCapacity << std::endl;
            UTIL_THROW("Invalid atom id");
         }
         file >> atomPtr->typeId;
         if (hasMolecules_) {
            file >> atomPtr->speciesId;
            if (atomPtr->speciesId < 0) {
               std::cout << "species Id  =" << atomPtr->speciesId << std::endl;
               UTIL_THROW("Negative species id");
            }
            file >> atomPtr->moleculeId; 
            if (atomPtr->moleculeId < 0) {
               std::cout << "molecule Id =" << atomPtr->moleculeId << std::endl;
               UTIL_THROW("Negative molecule id");
            }
            file >> atomPtr->atomId;
            if (atomPtr->atomId < 0) {
               std::cout << "atom id     =" << atomPtr->atomId << std::endl;
               UTIL_THROW("Negative atom id in molecule");
            }
         }
         file >> atomPtr->position;
         file >> atomPtr->velocity;

         // Finalize addition of new atom
         configuration().atoms().add();
      }

      // Read Covalent Groups
      #ifdef INTER_BOND
      if (configuration().bonds().capacity()) {
         readGroups(file, "BONDS", "nBond", configuration().bonds());
         //if (maskPolicy == MaskBonded) {
         //   setAtomMasks();
         //}
      }
      #endif

      // Optionally add atoms to species
      if (configuration().nSpecies() > 0) {
         if (!hasMolecules_) {
            UTIL_THROW("No atom context info in chosen ConfigIo");
         }
         int speciesId;
         AtomStorage::Iterator iter;
         configuration().atoms().begin(iter);
         for ( ; iter.notEnd(); ++iter) {
            speciesId = iter->speciesId;
            configuration().species(speciesId).addAtom(*iter);
         }
         for (int i = 0; configuration().nSpecies(); ++i) {
            configuration().species(i).isValid();
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
      file << configuration().boundary() << std::endl;
      file << std::endl;

      // Atoms
      int nAtom = configuration().atoms().size();
      file << "ATOMS" << std::endl;
      file << "nAtom" << Int(nAtom, 10) << std::endl;

      // Write atoms
      Vector r;
      AtomStorage::Iterator iter;
      configuration().atoms().begin(iter);
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
      if (configuration().bonds().capacity()) {
         writeGroups(file, "BONDS", "nBond", configuration().bonds());
      }
      #endif

   }
 
}
#endif
