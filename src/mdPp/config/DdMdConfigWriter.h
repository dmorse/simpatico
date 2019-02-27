#ifndef MDPP_DDMD_CONFIG_WRITER_H
#define MDPP_DDMD_CONFIG_WRITER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mdPp/config/ConfigWriter.h>
#include <mdPp/storage/GroupStorage.h>
#include <mdPp/storage/Configuration.h>

#include <util/format/Int.h>

namespace MdPp
{

   class Configuration;

   using namespace Util;

   /**
   * Native / default DdMd format for configuration files.
   *
   * \ingroup MdPp_ConfigWriter_Module
   */
   class DdMdConfigWriter  : public ConfigWriter
   {

   public:

      /**
      * Default constructor.
      *
      * \param writeAtomContexts Should atom contexts be written, if available.
      */
      DdMdConfigWriter(bool writeAtomContexts = true);

      /**
      * Constructor.
      *
      * \param configuration parent Configuration object.
      * \param writeAtomContexts true if file format has DdMd::AtomContext info
      */
      DdMdConfigWriter(Configuration& configuration, bool writeAtomContexts = false);

      /**
      * Write configuration file in DdMd default format.
      *
      * \param file output file stream (must be open on master)
      */
      virtual void writeConfig(std::ofstream& file);
   
   private:

      /*
      * Should atom context data be written, if available?
      */ 
      bool writeAtomContexts_;

      /**
      * Write BONDS, ANGLES, or DIHEDRALS data block.
      */
      template <int N>
      int writeGroups(std::ofstream& file, const char* sectionLabel, 
                      const char* nGroupLabel, GroupStorage<N>& groups);

   };

}
#endif
