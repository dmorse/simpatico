/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SmpConfigWriter.h"

#include <mdPp/chemistry/Atom.h>
#include <mdPp/chemistry/Group.h>
#include <mdPp/storage/SpeciesStorage.h>
#include <mdPp/storage/Configuration.h>

#include <util/misc/FlagSet.h>
#include <util/space/Vector.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

namespace MdPp
{

   using namespace Util;

   /*
   * Constructor.
   */
   SmpConfigWriter::SmpConfigWriter(Configuration& configuration)
    : ConfigWriter(configuration)
   {  setClassName("SmpConfigWriter"); }

   /* 
   * Write the configuration file.
   */
   void SmpConfigWriter::writeConfig(std::ofstream& file)
   {
      // Precondition
      if (!file.is_open()) {  
         UTIL_THROW("Error: File is not open"); 
      }

      using std::endl;

      // Write SPECIES block (structure of molecular species)
      int nSpecies = configuration().nSpecies();
      if (nSpecies > 0) {
         file << "SPECIES";
         file << endl << "nSpecies  " << nSpecies;
         int iSpecies, nMolecule, nAtom;
         int nAtomTot = 0;
         SpeciesStorage* speciesPtr;
         file << endl;
         for (iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
            file << endl << "species " << iSpecies;
            speciesPtr = &configuration().species(iSpecies);
            nMolecule = speciesPtr->size();
            file << endl << "  nMolecule  " << nMolecule << endl;
            speciesPtr->writeStructure(file, "  ");
            nAtom = speciesPtr->nAtom();
            nAtomTot += nMolecule*nAtom;
            file << endl;
         }
      }

      // Write BOUNDARY block (boundary dimensions)
      file << endl << "BOUNDARY";
      file << endl << configuration().boundary() << endl;
      file << std::endl;

      // Write ATOMS header
      file << endl << "ATOMS";
      file << endl << "format itmpv";
      int nAtom = configuration().atoms().size();
      file << "nAtom" << Int(nAtom, 10) << std::endl;

      // Write all atoms
      Vector r;
      AtomStorage::Iterator iter;
      configuration().atoms().begin(iter);
      for ( ; iter.notEnd(); ++iter) {
         file << endl;
         file << Int(iter->id, 10) 
              << Int(iter->typeId, 6);
         //if (hasMolecules_) {
         //   file << Int(iter->speciesId, 8) 
         //        << Int(iter->moleculeId, 8)
         //        << Int(iter->atomId, 8);
         //}
         r = iter->position;
         file << "\n" << r 
              << "\n" << iter->velocity << "\n";
      }

      // Write the groups
      #ifdef SIMP_BOND
      if (configuration().bonds().capacity()) {
         writeGroups(file, "BONDS", "nBond", configuration().bonds());
      }
      #endif
      #ifdef SIMP_ANGLE
      if (configuration().angles().capacity()) {
         writeGroups(file, "ANGLES", "nAngle", configuration().angles());
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (configuration().dihedrals().capacity()) {
         writeGroups(file, "DIHEDRALS", "nDihedral", 
                     configuration().dihedrals());
      }
      #endif

   }

   /*
   * Private method to write Group<N> objects.
   */
   template <int N>
   int SmpConfigWriter::writeGroups(
                  std::ofstream& file, 
                  const char* sectionLabel,
                  const char* nGroupLabel, 
                  GroupStorage<N>& groups)
   {
      ArrayIterator< Group<N> > iter;
      int nGroup = configuration().bonds().size();

      file << std::endl;
      file << sectionLabel << std::endl;
      file << nGroupLabel << Int(nGroup, 10) << std::endl;
      for (groups.begin(iter); iter.notEnd(); ++iter) {
         file << *iter << std::endl;
      }
      return nGroup;
   }

}
