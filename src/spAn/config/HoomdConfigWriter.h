#ifndef SPAN_HOOMD_CONFIG_WRITER_H
#define SPAN_HOOMD_CONFIG_WRITER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <spAn/config/ConfigWriter.h>
#include <spAn/storage/GroupStorage.h>
#include <spAn/storage/Configuration.h>

#include <util/format/Int.h>

namespace SpAn
{

   class Configuration;

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
   
   private:

      template <int N>
      int writeGroups(std::ofstream& file, const char* sectionLabel, 
                      const char* nGroupLabel, GroupStorage<N>& groups);

   };

   // Member functions templates

   /*
   * Private method to write Group<N> objects.
   */
   template <int N>
   int HoomdConfigWriter::writeGroups(std::ofstream& file, const char* sectionLabel,
                  const char* nGroupLabel, GroupStorage<N>& groups)
   {
      Configuration::BondIterator iter;
      int nGroup = configuration().bonds().size();

      file << std::endl;
      //file << sectionLabel << std::endl;
      //file << nGroupLabel << Int(nGroup, 10) << std::endl;
      for (groups.begin(iter); iter.notEnd(); ++iter) {
         file << *iter << std::endl;
      }
      return nGroup;
   }

}
#endif
