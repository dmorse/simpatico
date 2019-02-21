/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SmpConfigReader.h"

#include <mdPp/chemistry/Atom.h>
#include <mdPp/chemistry/Group.h>
#include <mdPp/storage/SpeciesStorage.h>
#include <mdPp/storage/Configuration.h>
#include <util/space/IntVector.h>
#include <util/param/Label.h>
#include <util/param/OptionalLabel.h>
#include <util/misc/FlagSet.h>

namespace MdPp
{

   using namespace Util;

   /*
   * Constructor.
   */
   SmpConfigReader::SmpConfigReader(Configuration& configuration)
    : ConfigReader(configuration)
   {  setClassName("SmpConfigReader"); }

   /* 
   * Read the configuration file.
   */
   void SmpConfigReader::readConfig(std::ifstream& file)
   {
      // Precondition
      if (!file.is_open()) {  
            UTIL_THROW("Error: File is not open"); 
      }
     
      // Make sure Label static buffer is clean on entry
      UTIL_CHECK(Label::isClear());

      // Read optional SPECIES block
      int nAtomTot = 0;
      OptionalLabel speciesLabel("SPECIES"); 
      bool hasSpecies = speciesLabel.match(file);
      if (hasSpecies) {

         // Read nSpecies, allocate array of species
         int nSpecies;
         file >> Label("nSpecies") >> nSpecies;
         UTIL_CHECK(nSpecies > 0);
         if (configuration().nSpecies() == 0) {
            configuration().setNSpecies(nSpecies);
         }
         UTIL_CHECK(nSpecies == configuration().nSpecies());

         // Loop over species 
         SpeciesStorage* speciesPtr;
         int iSpeciesIn, nMolecule;
         Label speciesLabel("species");
         Label nMoleculeLabel("nMolecule");
         for (int iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
            file >> Label("species") >> iSpeciesIn;
            UTIL_CHECK(iSpeciesIn == iSpecies);
            speciesPtr = &configuration().species(iSpecies);
            file >> Label("nMolecule") >> nMolecule;
            UTIL_CHECK(nMolecule >= 0);
            speciesPtr->setCapacity(nMolecule);
            speciesPtr->readStructure(file);
            speciesPtr->initialize();
            nAtomTot += nMolecule*speciesPtr->nAtom();
         }
      }

      // Read BOUNDARY block
      file >> Label("BOUNDARY");
      file >> configuration().boundary();

      // Read ATOMs block header
      file >> Label("ATOMS");
      Label orderedLabel("ordered", false); // optional label
      bool isOrdered = orderedLabel.match(file);
      std::string formatString;
      file >> Label("format") >> formatString;
      UTIL_CHECK(formatString.size() > 0);

      // Parse atom format string
      FlagSet atomFormat("itmpvs");
      atomFormat.setActualOrdered(formatString);
      bool hasAtomIndex = atomFormat.isActive('i');
      bool hasAtomTypeId = atomFormat.isActive('t');
      bool hasAtomContext = atomFormat.isActive('m');
      bool hasAtomPosition = atomFormat.isActive('p');
      bool hasAtomVelocity = atomFormat.isActive('v');
      bool hasAtomShift = atomFormat.isActive('s');

      // Consistency checks
      if (!isOrdered) {
         UTIL_CHECK(hasAtomIndex);
      }
      UTIL_CHECK(hasAtomPosition);
      if (hasSpecies) {
         UTIL_CHECK(isOrdered || hasAtomContext);
      }
      //if (hasAtomContext) {
      //   UTIL_CHECK(hasSpecies);
      //}
 
      // Read nAtom and allocate if necessary
      int nAtom;
      file >> Label("nAtom") >> nAtom;
      UTIL_CHECK(nAtom > 0);
      if (hasSpecies) {
         UTIL_CHECK(nAtom == nAtomTot);
      }
      if (configuration().atoms().capacity() == 0) {
         configuration().atoms().allocate(nAtom);
      } 
      UTIL_CHECK(configuration().atoms().capacity() >= nAtom);

      // Read atoms 
      Atom* atomPtr = 0;
      IntVector shift;
      int atomCapacity = configuration().atoms().capacity(); 
      for (int i = 0; i < nAtom; ++i) {

         // Get pointer to new atom 
         atomPtr = configuration().atoms().newPtr();
 
         if (hasAtomIndex) {
            file >> atomPtr->id;
            UTIL_CHECK(atomPtr->id >= 0);
            if (isOrdered) {
               UTIL_CHECK(atomPtr->id == i);
            }
            UTIL_CHECK(atomPtr->id < atomCapacity);
         }
         if (hasAtomTypeId) {
            file >> atomPtr->typeId;
            UTIL_CHECK(atomPtr->typeId >= 0);
         }
         if (hasAtomContext) {
            file >> atomPtr->speciesId;
            UTIL_CHECK(atomPtr->speciesId >= 0);
            file >> atomPtr->moleculeId;
            UTIL_CHECK(atomPtr->moleculeId >= 0);
            file >> atomPtr->atomId;
            UTIL_CHECK(atomPtr->atomId >= 0);
         } 
         file >> atomPtr->position;
         if (hasAtomVelocity) {
            file >> atomPtr->velocity;
         }
         if (hasAtomShift) {
            file >> shift;
            #ifdef MDPP_SHIFT
            atomPtr->shift() = shift;
            #endif
         }

         // Shift atom positions to primary cell
         #ifdef MDPP_SHIFT
         configuration().boundary().shift(atomPtr->position, 
                                          atomPtr->shift);
         #else
         configuration().boundary().shift(atomPtr->position);
         #endif

         // Finalize addition of new atom
         configuration().atoms().add();
      }

      // If hasSpecies, add atoms to species
      if (hasSpecies) {
         if (hasAtomContext) {
            addAtomsToSpecies();
         } else {
            bool success;
            success = setAtomContexts();
            if (success) {
               addAtomsToSpecies();
            }
         }
         // TODO: Make groups based on species description
      } else {
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
         #ifdef SIMP_IMPROPERS
         readGroups(file, "IMPROPERS", "nDihedral", 
                    configuration().dihedrals());
         #endif
      }

      // Make sure static Label buffer is clean on exit
      UTIL_CHECK(Label::isClear());
   }

   /*
   * Private method to read Group<N> objects.
   */
   template <int N>
   int SmpConfigReader::readGroups(std::ifstream& file, 
                  const char* sectionLabel,
                  const char* nGroupLabel,
                  GroupStorage<N>& groups)
   {
      // Check for sectionLabel heading
      OptionalLabel speciesLabel(sectionLabel); 
      bool hasGroups = speciesLabel.match(file);

      int nGroup = 0;  
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
