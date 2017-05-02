/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "HoomdConfigWriter.h"

#include <tools/chemistry/Atom.h>
#include <tools/chemistry/Group.h>
#include <tools/chemistry/Species.h>
//#include <tools/chemistry/MaskPolicy.h>
#include <tools/storage/GroupStorage.h>
#include <tools/storage/Configuration.h>

#include <util/space/Vector.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>
#include <util/misc/ioUtil.h>

namespace Tools
{

   using namespace Util;

   /*
   * Constructor.
   */
   HoomdConfigWriter::HoomdConfigWriter(Configuration& configuration)
    : ConfigWriter(configuration, true),
      hasTypeMaps_(false)
   {  setClassName("HoomdConfigWriter"); }

   /*
   * Read auxiliary type map file.
   */
   void HoomdConfigWriter::readAuxiliaryFile(std::ifstream& file)
   {
      bool notEnd;
      std::stringstream line;

      notEnd = getNextLine(file, line);
      if (notEnd) {
         checkString(line, "ATOM");
         checkString(line, "TYPES:");
         atomTypeMap_.read(file);
      }

      notEnd = getNextLine(file, line);
      if (notEnd) {
         checkString(line, "BOND");
         checkString(line, "TYPES:");
         bondTypeMap_.read(file);
      }

      notEnd = getNextLine(file, line);
      if (notEnd) {
         checkString(line, "ANGLE");
         checkString(line, "TYPES:");
         angleTypeMap_.read(file);
      }

      hasTypeMaps_ = true; 
   }

   /*
   * Private function template for writing groups.
   */
   template <int N> void
   HoomdConfigWriter::writeGroups(std::ofstream& file, 
                                  const std::string& label,
                                  const GroupStorage<N>& storage,
                                  const TypeMap& map)
   {
      int nGroup = storage.size();
      file << "<" << label 
           << " num=\"" << nGroup << "\">\n";
      std::string typeName;
      int j;
      ConstArrayIterator< Group<N> > iter;
      for (storage.begin(iter); iter.notEnd(); ++iter) {
         typeName = map.name(iter->typeId);
         file << typeName << " ";
         for (j = 0; j < N; ++j) {
            file << iter->atomIds[j] << " ";
         }            
         file << "\n";
      }
      file << "</" << label << ">\n";
   }

   /* 
   * Write the configuration file.
   */
   void HoomdConfigWriter::writeConfig(std::ofstream& file)
   {
      // Precondition
      if (!file.is_open()) {  
         UTIL_THROW("Error: File is not open"); 
      }
      if (!hasTypeMaps_) {  
         UTIL_THROW("Error: Must read type map auxiliary file before config"); 
      }
      int nAtom = configuration().atoms().size();

      file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
      file << "<hoomd_xml version=\"1.5\">\n";
      file << "<configuration time_step=\"0\" dimensions=\"3\" natoms=\""
           << nAtom << "\" >\n";

      // Set ostream floating point precision, restore default format
      file.precision(12);
      file.unsetf(std::ios_base::floatfield);

      // Write box
      Boundary& boundary = configuration().boundary();
      Vector lengths = boundary.lengths();
      file << "<box"
           << " lx=\"" << lengths[0] << "\""
           << " ly=\"" << lengths[1] << "\""
           << " lz=\"" << lengths[2] << "\""
           << " xy=\"0\""
           << " xz=\"0\""
           << " yz=\"0\""
           << "/>\n";

      // Write atom position node
      file << "<position num=\"" << nAtom << "\">\n";
      Vector rg, rc;
      AtomStorage::Iterator iter;
      int j;
      configuration().atoms().begin(iter);
      for (; iter.notEnd(); ++iter) {
         boundary.shift(iter->position);
         boundary.transformCartToGen(iter->position, rg);
         for (j = 0; j < Dimension; ++j) {
            if (rg[j] >= 0.5) {
               rg[j] -= 1.0;
            }
         }
         boundary.transformGenToCart(rg, rc);
         for (j = 0; j < Dimension; ++j) {
            file.precision(12);
            file << rc[j] << " ";
         }
         file << "\n";
      }
      file << "</position>\n";

      #if 0
      // Write velocity
      file << "<velocity num=\"" << nAtom << "\">\n";
      configuration().atoms().begin(iter);
      for (; iter.notEnd(); ++iter) {
         file << iter->velocity << "\n";
      }
      file << "</velocity>\n";
      #endif

      // Write type
      file << "<type num=\"" << nAtom << "\">\n";
      configuration().atoms().begin(iter);
      for (; iter.notEnd(); ++iter) {
         file << atomTypeMap_.name(iter->typeId) << "\n";
      }
      file << "</type>\n";

      // Write covalent groups, as needed
      #ifdef SIMP_BOND
      if (configuration().bonds().size()) {
         writeGroups(file, "bond", configuration().bonds(), 
                     bondTypeMap_);
      }
      #endif
      #ifdef SIMP_ANGLE
      if (configuration().angles().size()) {
         writeGroups(file, "angle", configuration().angles(), 
                     angleTypeMap_);
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (configuration().dihedrals().size()) {
         writeGroups(file, "dihedral", configuration().dihedrals(), 
                     dihedralTypeMap_);
      }
      #endif

      file << "</configuration>\n";
      file << "</hoomd_xml>\n";
   }
 
}
