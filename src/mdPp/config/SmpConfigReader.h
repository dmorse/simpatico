#ifndef MDPP_SMP_CONFIG_READER_H
#define MDPP_SMP_CONFIG_READER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mdPp/config/ConfigReader.h>
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
   * \ingroup MdPp_ConfigReader_Module
   */
   class SmpConfigReader  : public ConfigReader
   {

   public:

      /**
      * Default constructor.
      */
      SmpConfigReader();

      /**
      * Constructor.
      *
      * \param configuration parent Configuration object.
      */
      SmpConfigReader(Configuration& configuration);

      /**
      * Read configuration file in DdMd default format.
      *
      * \param file input file stream
      */
      virtual void readConfig(std::ifstream& file);

   private:

      bool hasMolecules_;

      template <int N>
      int readGroups(std::ifstream& file, const char* sectionLabel, 
                     const char* nGroupLabel, GroupStorage<N>& groups);

   };

}
#endif
