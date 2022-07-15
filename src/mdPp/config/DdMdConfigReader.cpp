/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "DdMdConfigReader.h"

#include <mdPp/chemistry/Atom.h>
#include <mdPp/chemistry/Group.h>
#include <mdPp/storage/Configuration.h>
#include <util/param/OptionalLabel.h>

namespace MdPp
{

   using namespace Util;

   /*
   * Constructor.
   */
   DdMdConfigReader::DdMdConfigReader(Configuration& configuration, 
                                      bool hasMolecules)
    : ConfigReader(configuration),
      hasAtomContexts_(hasMolecules)
   {  setClassName("DdMdConfigReader"); }

   /*
   * Read a configuration file.
   */
   void DdMdConfigReader::readConfig(std::ifstream& file)
   {
      // Precondition
      if (!file.is_open()) {  
            UTIL_THROW("Error: File is not open"); 
      }

      using std::endl;

      // Read and broadcast boundary
      file >> Label("BOUNDARY");
      file >> configuration().boundary();

      // Read ATOMS header, allocate if necessary
      file >> Label("ATOMS");
      int nAtom;          
      file >> Label("nAtom") >> nAtom;
      UTIL_CHECK(nAtom > 0);
      if (configuration().atoms().capacity() == 0) {
         configuration().atoms().allocate(nAtom); 
      }
      int atomCapacity = configuration().atoms().capacity(); 
      UTIL_CHECK(nAtom <= atomCapacity);
      configuration().setHasAtomContexts(hasAtomContexts_);

      // Read and distribute atoms
      Atom* atomPtr;
      for (int i = 0; i < nAtom; ++i) {

         // Get pointer to new atom 
         atomPtr = configuration().atoms().newPtr();
 
         file >> atomPtr->id;
         if (atomPtr->id < 0) {
            std::cout << "atom id =" << atomPtr->id << endl;
            UTIL_THROW("Negative atom id");
         }
         if (atomPtr->id >= atomCapacity) {
            std::cout << "atom id      =" << atomPtr->id << endl;
            std::cout << "atomCapacity =" << atomCapacity << endl;
            UTIL_THROW("Invalid atom id");
         }
         file >> atomPtr->typeId;
         if (hasAtomContexts_) {
            file >> atomPtr->speciesId;
            if (atomPtr->speciesId < 0) {
               std::cout << "species Id  =" << atomPtr->speciesId << endl;
               UTIL_THROW("Negative species id");
            }
            file >> atomPtr->moleculeId; 
            if (atomPtr->moleculeId < 0) {
               std::cout << "molecule Id =" << atomPtr->moleculeId << endl;
               UTIL_THROW("Negative molecule id");
            }
            file >> atomPtr->atomId;
            if (atomPtr->atomId < 0) {
               std::cout << "atom id     =" << atomPtr->atomId << endl;
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
      readGroups(file, "BONDS", "nBond", configuration().bonds());
      #endif

      #ifdef SIMP_ANGLE
      readGroups(file, "ANGLES", "nAngle", configuration().angles());
      #endif

      #ifdef SIMP_DIHEDRAL
      readGroups(file, "DIHEDRALS", "nDihedral", 
                 configuration().dihedrals());
      #endif

   }
 
   /*
   * Private function template to read Group<N> objects.
   */
   template <int N>
   int DdMdConfigReader::readGroups(std::ifstream& file, 
                  const char* sectionLabel,
                  const char* nGroupLabel,
                  GroupStorage<N>& groups)
   {
      int nGroup = 0;  

      // Check for sectionLabel heading
      OptionalLabel speciesLabel(sectionLabel); 
      bool hasGroups = speciesLabel.match(file);

      if (hasGroups) {

         // Read nGroup, allocate if necessary
         file >> Label(nGroupLabel) >> nGroup;
         UTIL_CHECK(nGroup > 0);
         if (groups.capacity() == 0) {
            groups.allocate(nGroup);
         }
         UTIL_CHECK(groups.capacity() >= nGroup);

         // Read groups
         Group<N>* groupPtr;
         for (int i = 0; i < nGroup; ++i) {
            groupPtr = groups.newPtr();
            file >> *groupPtr;
         }

      }
      return nGroup;
   }

}
