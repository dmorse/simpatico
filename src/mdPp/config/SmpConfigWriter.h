#ifndef MDPP_SMP_CONFIG_WRITER_H
#define MDPP_SMP_CONFIG_WRITER_H

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
   * Simpatico configuration file writer.
   *
   * \ingroup MdPp_ConfigWriter_Module
   */
   class SmpConfigWriter  : public ConfigWriter
   {

   public:

      /**
      * Default constructor.
      */
      SmpConfigWriter();

      /**
      * Constructor.
      *
      * \param configuration parent Configuration object.
      */
      SmpConfigWriter(Configuration& configuration);

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

}
#endif
