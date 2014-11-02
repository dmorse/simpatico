#ifndef TOOLS_HOOMD_CONFIG_WRITER_H
#define TOOLS_HOOMD_CONFIG_WRITER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <tools/config/ConfigWriter.h>       // base class
#include <tools/config/TypeMap.h>            // member

namespace Tools
{

   class Configuration;
   template <int N> class GroupStorage;

   using namespace Util;

   /**
   * Hoomd-blue XML format for configuration files.
   *
   * \ingroup Tools_ConfigWriter_Module
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
      TypeMap improperTypeMap_;
      bool hasTypeMaps_;

      template <int N>
      void writeGroups(std::ofstream& file,
                       const std::string& label,
                       const GroupStorage<N>& groups,
                       const TypeMap& map);

   };

}
#endif
