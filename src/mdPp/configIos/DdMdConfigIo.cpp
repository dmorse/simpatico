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
#include <mdPp/processor/Processor.h>
//#include <mdPp/chemistry/MaskPolicy.h>

#include <util/space/Vector.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

namespace MdPp
{

   using namespace Util;

   /*
   * Constructor.
   */
   DdMdConfigIo::DdMdConfigIo(Processor& processor, bool hasMolecules)
    : ConfigIo(processor),
      hasMolecules_(hasMolecules)
   {  setClassName("DdMdConfigIo"); }

   /*
   * Private method to read Group<N> objects.
   */
   int DdMdConfigIo::readBonds(std::ifstream& file, 
                  const char* sectionLabel,
                  const char* nGroupLabel)
   {
      int nGroup;  // Total number of groups in file
      file >> Label(sectionLabel);
      file >> Label(nGroupLabel) >> nGroup;
      Group<2>* groupPtr;
      for (int i = 0; i < nGroup; ++i) {
         groupPtr = processor().newBondPtr();
         file >> *groupPtr;
      }
      return nGroup;
   }

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
      file >> processor().boundary();

      // Read and distribute atoms

      // Read atoms
      Atom* atomPtr;
      int atomCapacity = processor().atomCapacity(); // Maximum allowed id + 1
      int nAtom;          
      file >> Label("ATOMS");
      file >> Label("nAtom") >> nAtom;
      for (int i = 0; i < nAtom; ++i) {

         // Get pointer to new atom 
         atomPtr = processor().newAtomPtr();
 
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

         processor().addAtom();
      }

      // Read Covalent Groups
      #ifdef INTER_BOND
      if (processor().bondCapacity()) {
         readBonds(file, "BONDS", "nBond");
         //if (maskPolicy == MaskBonded) {
         //   setAtomMasks();
         //}
      }
      #endif
   }

   int DdMdConfigIo::writeBonds(std::ofstream& file, 
                  const char* sectionLabel,
                  const char* nGroupLabel)
   {
      Processor::BondIterator iter;
      int nGroup = processor().nBond();

      file << std::endl;
      file << sectionLabel << std::endl;
      file << nGroupLabel << Int(nGroup, 10) << std::endl;
      for (processor().initBondIterator(iter); iter.notEnd(); ++iter) {
         file << *iter << std::endl;
      }
      return nGroup;
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
      file << processor().boundary() << std::endl;
      file << std::endl;

      // Atoms
      int nAtom = processor().nAtom();
      file << "ATOMS" << std::endl;
      file << "nAtom" << Int(nAtom, 10) << std::endl;

      // Write atoms
      Vector r;
      Processor::AtomIterator iter;
      processor().initAtomIterator(iter);
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
      if (processor().bondCapacity()) {
         writeBonds(file, "BONDS", "nBond");
      }
      #endif

   }
 
}
#endif
