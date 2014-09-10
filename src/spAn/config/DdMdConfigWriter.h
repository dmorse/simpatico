#ifndef SPAN_DDMD_CONFIG_WRITER_H
#define SPAN_DDMD_CONFIG_WRITER_H

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
   * Native / default DdMd format for configuration files.
   *
   * \ingroup SpAn_ConfigWriter_Module
   */
   class DdMdConfigWriter  : public ConfigWriter
   {

   public:

      /**
      * Default constructor.
      */
      DdMdConfigWriter(bool hasMolecules = false);

      /**
      * Constructor.
      *
      * \param configuration parent Configuration object.
      * \param hasMolecules true if file format has DdMd::AtomContext info
      */
      DdMdConfigWriter(Configuration& configuration, bool hasMolecules = false);

      /**
      * Write configuration file in DdMd default format.
      *
      * \param file output file stream (must be open on master)
      */
      virtual void writeConfig(std::ofstream& file);
   
   private:

      bool hasMolecules_;

      template <int N>
      int writeGroups(std::ofstream& file, const char* sectionLabel, 
                      const char* nGroupLabel, GroupStorage<N>& groups);

   };

   // Member functions templates

   /*
   * Private method to write Group<N> objects.
   */
   template <int N>
   int DdMdConfigWriter::writeGroups(std::ofstream& file, const char* sectionLabel,
                  const char* nGroupLabel, GroupStorage<N>& groups)
   {
      Configuration::BondIterator iter;
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
#endif
