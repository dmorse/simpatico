#ifndef DDMD_DDMD_CONFIG_IO_H
#define DDMD_DDMD_CONFIG_IO_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/configIos/ConfigIo.h>

namespace DdMd
{

   class Simulation;

   using namespace Util;

   /**
   * Native / default format for configuration files.
   *
   * DdMdConfigIo is a ConfigIo that implements the default
   * configuration file format for ddSim.
   *
   * \ingroup DdMd_ConfigIo_Module
   */
   class DdMdConfigIo  : public ConfigIo
   {

   public:

      /**
      * Default constructor.
      */
      DdMdConfigIo();

      /**
      * Constructor.
      *
      * \param simulation parent Simulation object.
      */
      DdMdConfigIo(Simulation& simulation);

      /**
      * Read configuration file.
      *
      * This routine opens and reads a file on the master, and distributes
      * atom data among the processors.
      *
      * \param file input file stream
      * \param maskPolicy MaskPolicy to be used in setting atom masks
      */
      virtual void readConfig(std::ifstream& file, MaskPolicy maskPolicy);

      /**
      * Write configuration file.
      *
      * This routine opens and writes a file on the master,
      * collecting atom data from all processors.
      *
      * \param file output file stream
      */
      virtual void writeConfig(std::ofstream& file);
   
   private:

      /**
      * Read Group<N> objects from file. 
      */
      template <int N>
      int readGroups(std::ifstream& file, 
                     const char* sectionLabel, const char* nGroupLabel,
                     GroupDistributor<N>& distributor);

      /**
      * Write Group<N> objects to file. 
      */
      template <int N>
      int writeGroups(std::ofstream& file, 
                      const char* sectionLabel, const char* nGroupLabel,
                      GroupStorage<N>& storage, GroupCollector<N>& collector);
   
   };

}
#endif
