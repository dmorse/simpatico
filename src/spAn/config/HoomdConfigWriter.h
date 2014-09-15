#ifndef SPAN_HOOMD_CONFIG_WRITER_H
#define SPAN_HOOMD_CONFIG_WRITER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <spAn/config/ConfigWriter.h>       // base class
#include <span/config/TypeMap.h>            // member

namespace SpAn
{

   class Configuration;
   template <int N> class GroupStorage;

   using namespace Util;

   /**
   * Hoomd-blue XML format for configuration files.
   *
   * \ingroup SpAn_ConfigWriter_Module
   */
   class HoomdConfigWriter  : public ConfigWriter
   {

   public:

      /**
      * Default constructor.
      */
      HoomdConfigWriter();

      /**
      * Constructor.
      *
      * \param configuration parent Configuration object.
      */
      HoomdConfigWriter(Configuration& configuration);

      /**
      * Write configuration file in DdMd default format.
      *
      * \param file output file stream (must be open on master)
      */
      virtual void writeConfig(std::ofstream& file);
   
      /**
      * Read auxiliary file with type map information.
      *
      * \param file input file stream
      */
      virtual void readAuxiliaryFile(std::ifstream& file);

   private:

      TypeMap atomTypeMap_;
      TypeMap bondTypeMap_;
      TypeMap angleTypeMap_;
      TypeMap dihedralTypeMap_;

      template <int N>
      void writeGroups(std::ofstream& file,
                       const std::string& label,
                       const GroupStorage<N>& groups,
                       const TypeMap& map);

   };

   #if 0
   // Member functions templates

   /*
   * Private method to write Group<N> objects.
   */
   template <int N> void
   HoomdConfigWriter::writeGroups(std::ofstream& file, 
                                  const std::string& label,
                                  const GroupStorage<N>& groups,
                                  const TypeMap& map)
   {
      file << "<" << label << ">\n";
      std::string typeName;
      int j;
      ConstArrayIterator< Group<N> > iter;
      for (groups.begin(iter); iter.notEnd(); ++iter) {
         typeName = map.name(iter->typeId);
         file << typeName << " ";
         for (j = 0; j < N; ++j) {
            file << iter->atomIds[j] << " ";
         }            
         file << "\n";
      }
      file << "</" << label << ">\n\n";
   }
   #endif

}
#endif
