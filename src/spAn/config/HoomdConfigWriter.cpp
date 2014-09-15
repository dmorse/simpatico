#ifndef SPAN_HOOMD_CONFIG_WRITER_CPP
#define SPAN_HOOMD_CONFIG_WRITER_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "HoomdConfigWriter.h"

#include <spAn/chemistry/Atom.h>
#include <spAn/chemistry/Group.h>
#include <spAn/chemistry/Species.h>
//#include <spAn/chemistry/MaskPolicy.h>
#include <spAn/storage/GroupStorage.h>
#include <spAn/storage/Configuration.h>

#include <util/space/Vector.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>
#include <util/misc/ioUtil.h>

namespace SpAn
{

   using namespace Util;

   /*
   * Constructor.
   */
   HoomdConfigWriter::HoomdConfigWriter(Configuration& configuration)
    : ConfigWriter(configuration, true)
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
      file << "</" << label << ">\n\n";
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

      file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
      file << "<hoomd_xml version=\"1.5\">\n\n";
      file << "<configuration time_step=\"0\">\n\n";

      // Write box
      Vector lengths = configuration().boundary().lengths();
      file << "<box"
           << " lx=\"" << lengths[0] << "\""
           << " ly=\"" << lengths[1] << "\""
           << " lz=\"" << lengths[1] << "\""
           << " />\n\n";

      // Write position
      int nAtom = configuration().atoms().size();
      file << "<position num=\"" << nAtom << "\">\n";
      AtomStorage::Iterator iter;
      configuration().atoms().begin(iter);
      for (; iter.notEnd(); ++iter) {
         file << iter->position << "\n";
      }
      file << "</position>\n\n";

      // Write position
      file << "<velocity num=\"" << nAtom << "\">\n";
      configuration().atoms().begin(iter);
      for (; iter.notEnd(); ++iter) {
         file << iter->position << "\n";
      }
      file << "</velocity>\n\n";

      // Write the groups
      #ifdef INTER_BOND
      if (configuration().bonds().size()) {
         writeGroups(file, "bond", configuration().bonds(), 
                     bondTypeMap_);
      }
      #endif

      #ifdef INTER_ANGLE
      if (configuration().angles().size()) {
         writeGroups(file, "angle", configuration().angles(), 
                     angleTypeMap_);
      }
      #endif

      #ifdef INTER_DIHEDRAL
      if (configuration().dihedrals().size()) {
         writeGroups(file, "dihedral", configuration().dihedrals(), 
                     dihedralTypeMap_);
      }
      #endif

      file << "</configuration>\n";
      file << "</hoomd_xml>\n";
   }
 
}
#endif
