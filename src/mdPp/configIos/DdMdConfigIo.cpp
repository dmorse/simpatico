#ifndef MDPP_DDMD_CONFIG_IO_CPP
#define MDPP_DDMD_CONFIG_IO_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "DdMdConfigIo.h"

#include <mdPp/Processor.h>
#include <mdPp/chemistry/Atom.h>
#include <mdPp/chemistry/Group.h>
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
   DdMdConfigIo::DdMdConfigIo(bool hasMolecules)
    : ConfigIo(),
      hasMolecules_(hasMolecules)
   {  setClassName("DdMdConfigIo"); }

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
         groupPtr = processor().addBond();
         file >> *groupPtr;
      }
      return nGroup;
   }

   /*
   * Read a configuration file.
   */
   void DdMdConfigIo::readConfig(std::ifstream& file)
   // , MaskPolicy maskPolicy)
   {
      // Precondition
      if (!file.is_open()) {  
            UTIL_THROW("Error: File is not open"); 
      }

      // Read and broadcast boundary
      file >> Label("BOUNDARY");
      file >> processor().boundary();

      // Atoms 
      int nAtom;  // Total number of atoms in file

      // Read and distribute atoms
      file >> Label("ATOMS");
      file >> Label("nAtom") >> nAtom;

      // int atomCapacity = processor().atomCapacity();

      // Read atoms
      Vector r;
      Atom* atomPtr;
      int id;
      int typeId;
      int aId;
      int mId;
      int sId;

      for (int i = 0; i < nAtom; ++i) {

         // Get pointer to new atom 
         atomPtr = processor().newAtomPtr();
 
         file >> id >> typeId;
         //if (id < 0 || id >= atomCapacity) {
         //   UTIL_THROW("Invalid atom id");
         //}
         atomPtr->id = id;
         atomPtr->typeId = typeId;
         if (hasMolecules_) {
            file >> sId >> mId >> aId;
            if (aId < 0) {
               UTIL_THROW("Invalid Atom");
            }
            if (mId < 0) {
               UTIL_THROW("Invalid Molecule");
            }
            if (sId < 0) {
               UTIL_THROW("Invalid Specie");
            }
            atomPtr->atomId = aId;
            atomPtr->moleculeId = mId;
            atomPtr->speciesId = sId;
         }
         file >> r;
         file >> atomPtr->velocity;

         processor().addAtom();
      }

      // Read Covalent Groups
      #if 0
      bool hasGhosts = false;
      #ifdef INTER_BOND
      if (processor().bondCapacity()) {
         readBonds(file, "BONDS", "nBond");
         // processor().isValid(atomStorage(), domain().communicator(), hasGhosts);
         // Set atom "mask" values
         //if (maskPolicy == MaskBonded) {
         //   setAtomMasks();
         //}
      }
      #endif
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

      #if 0
      // Write the groups
      #ifdef INTER_BOND
      if (processor().bondCapacity()) {
         writeBonds(file, "BONDS", "nBond", processor());
      }
      #endif
      #endif

   }
 
}
#endif
