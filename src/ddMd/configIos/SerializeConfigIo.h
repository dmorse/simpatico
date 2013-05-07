#ifndef DDMD_DDMD_SERIALIZE_CONFIG_IO_H
#define DDMD_DDMD_SERIALIZE_CONFIG_IO_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/configIos/ConfigIo.h>
#include <util/archives/Serializable.h>

namespace DdMd
{

   class Simulation;

   using namespace Util;

   /**
   * Save / load configuration to archives.
   *
   * \ingroup DdMd_ConfigIo_Module
   */
   class SerializeConfigIo  : public ConfigIo
   {

   public:

      /**
      * Default constructor.
      */
      SerializeConfigIo();

      /**
      * Constructor.
      *
      * \param simulation parent Simulation object.
      */
      SerializeConfigIo(Simulation& simulation);

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
   
      /**
      * Load configuration.
      *
      * \param ar input archive
      * \param maskPolicy MaskPolicy to be used in setting atom masks
      */
      void loadConfig(Serializable::IArchive& ar, MaskPolicy maskPolicy);

      /**
      * Save configuration.
      *
      * \param ar output archive
      */
      void saveConfig(Serializable::OArchive& ar);
   
   private:

      /**
      * Read Group<N> objects from file. 
      */
      template <int N>
      int loadGroups(Serializable::IArchive& ar, 
                     GroupDistributor<N>& distributor);

      /**
      * Write Group<N> objects to file. 
      */
      template <int N>
      int saveGroups(Serializable::OArchive& ar, 
                     GroupStorage<N>& storage, GroupCollector<N>& collector);
   
   };

}
#endif
