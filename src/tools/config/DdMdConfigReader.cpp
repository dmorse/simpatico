/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "DdMdConfigReader.h"

#include <tools/chemistry/Atom.h>
#include <tools/chemistry/Group.h>
#include <tools/chemistry/Species.h>
//#include <tools/chemistry/MaskPolicy.h>
#include <tools/storage/Configuration.h>

#include <util/space/Vector.h>
//#include <util/format/Int.h>
//#include <util/format/Dbl.h>

namespace Tools
{

   using namespace Util;

   /*
   * Constructor.
   */
   DdMdConfigReader::DdMdConfigReader(Configuration& configuration, bool hasMolecules)
    : ConfigReader(configuration),
      hasMolecules_(hasMolecules)
   {  setClassName("DdMdConfigReader"); }

   /*
   * Private method to read Group<N> objects.
   */
   template <int N>
   int DdMdConfigReader::readGroups(std::ifstream& file, 
                  const char* sectionLabel,
                  const char* nGroupLabel,
                  GroupStorage<N>& groups)
   {
      int nGroup;  // Total number of groups in file
      file >> Label(sectionLabel);
      file >> Label(nGroupLabel) >> nGroup;
      Group<N>* groupPtr;
      for (int i = 0; i < nGroup; ++i) {
         groupPtr = groups.newPtr();
         file >> *groupPtr;
      }
      return nGroup;
   }

   /*
   * Read a configuration file.
   */
   void DdMdConfigReader::readConfig(std::ifstream& file)
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
      #ifdef SIMP_BOND
      if (configuration().bonds().capacity()) {
         readGroups(file, "BONDS", "nBond", configuration().bonds());
         //if (maskPolicy == MaskBonded) {
         //   setAtomMasks();
         //}
      }
      #endif

      #ifdef SIMP_ANGLE
      if (configuration().angles().capacity()) {
         readGroups(file, "ANGLES", "nAngle", configuration().angles());
      }
      #endif

      #ifdef SIMP_DIHEDRAL
      if (configuration().dihedrals().capacity()) {
         readGroups(file, "DIHEDRALS", "nDihedral", configuration().dihedrals());
      }
      #endif

      // Optionally add atoms to species
      if (configuration().nSpecies() > 0) {
         if (hasMolecules_) {
            addAtomsToSpecies();
         } else {
            bool success;
            success = setAtomContexts();
            if (success) {
               addAtomsToSpecies();
            }
         }
      }

   }
 
}
