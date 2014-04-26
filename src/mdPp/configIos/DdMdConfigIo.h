#ifndef MDPP_DDMD_CONFIG_IO_H
#define MDPP_DDMD_CONFIG_IO_H

/*
* Simpatico - Processor Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mdPp/configIos/ConfigIo.h>

namespace MdPp
{

   class Processor;

   using namespace Util;

   /**
   * Native / default DdMd format for configuration files.
   *
   * \ingroup DdMd_ConfigIo_Module
   */
   class DdMdConfigIo  : public ConfigIo
   {

   public:

      /**
      * Default constructor.
      */
      DdMdConfigIo(bool hasMolecules = false);

      /**
      * Constructor.
      *
      * \param processor parent Processor object.
      */
      DdMdConfigIo(Processor& processor, bool hasMolecules = false);

      /**
      * Read configuration file.
      *
      * This routine opens and reads a file on the master, and distributes
      * atom data among the processors.
      *
      * \param file input file stream
      * \param maskPolicy MaskPolicy to be used in setting atom masks
      */
      virtual void readConfig(std::ifstream& file);

      /**
      * Write configuration file.
      *
      * This routine writes a file in DdMd default format.
      *
      * \param file output file stream (must be open on master)
      */
      virtual void writeConfig(std::ofstream& file);
   
   private:

      bool hasMolecules_;

      int readBonds(std::ifstream& file, 
                    const char* sectionLabel, const char* nGroupLabel);
      int writeBonds(std::ofstream& file, 
                    const char* sectionLabel, const char* nGroupLabel);

   };

}
#endif
