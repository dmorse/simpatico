/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "DdMdConfigWriter.h"

#include <tools/chemistry/Atom.h>
#include <tools/chemistry/Group.h>
#include <tools/chemistry/Species.h>
//#include <tools/chemistry/MaskPolicy.h>
#include <tools/storage/Configuration.h>

#include <util/space/Vector.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

namespace Tools
{

   using namespace Util;

   /*
   * Constructor.
   */
   DdMdConfigWriter::DdMdConfigWriter(Configuration& configuration, bool hasMolecules)
    : ConfigWriter(configuration),
      hasMolecules_(hasMolecules)
   {  setClassName("DdMdConfigWriter"); }

   /*
   * Private method to write Group<N> objects.
   */
   template <int N>
   int DdMdConfigWriter::writeGroups(std::ofstream& file, const char* sectionLabel,
                  const char* nGroupLabel, GroupStorage<N>& groups)
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

   /* 
   * Write the configuration file.
   */
   void DdMdConfigWriter::writeConfig(std::ofstream& file)
   {
      // Precondition
      if (!file.is_open()) {  
         UTIL_THROW("Error: File is not open"); 
      }

      // Write boundary box dimensions
      file << "BOUNDARY" << std::endl << std::endl;
      file << configuration().boundary() << std::endl;
      file << std::endl;

      // Write atoms
      file << "ATOMS" << std::endl;
      int nAtom = configuration().atoms().size();
      file << "nAtom" << Int(nAtom, 10) << std::endl;
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
         writeGroups(file, "DIHEDRALS", "nDihedral", configuration().dihedrals());
      }
      #endif

   }
 
}
