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
   * Save / load configuration from / to an archive.
   *
   * The virtual functions readConfig and writeConfig, which take file 
   * objects parameters, read from or write to files that stores archives 
   * of type Serializable::IArchive or Serializable::OArchive. The 
   * loadConfig and saveConfig methods take archive object arguments.
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
      * Call on all processors, but reads from file only on master.
      *
      * This function opens and reads a file on the master, and distributes
      * atom data among the processors. This calls loadConfig() internally,
      * and thus has all the preconditions of loadConfig().
      *
      * \param file input file stream (must be open on master)
      * \param maskPolicy MaskPolicy to be used in setting atom masks
      */
      virtual void readConfig(std::ifstream& file, MaskPolicy maskPolicy);

      /**
      * Write configuration file.
      *
      * Call on all processors, but file must be open only on master.
      *
      * This function opens and writes a file on the master, collecting
      * atom data from all processors. It calls saveConfig() internally.
      *
      * \param file output file stream (must be open on master)
      */
      virtual void writeConfig(std::ofstream& file);
   
      /**
      * Load configuration from an archive.
      *
      * Call on all processors, but save from archive only on master.
      *
      * \pre  There are no atoms, ghosts, or groups.
      * \pre  AtomStorage is set for scaled / generalized coordinates
      *
      * \post atomic coordinates are scaled / generalized
      * \post there are no ghosts
      *
      * \param ar input archive
      * \param maskPolicy MaskPolicy to be used in setting atom masks
      */
      void loadConfig(Serializable::IArchive& ar, MaskPolicy maskPolicy);

      /**
      * Save configuration.
      *
      * Call on all processors, but uses archive only on master.
      *
      * This function allows the AtomStorage to be set for either
      * generalized or Cartesian coordinates, but always saves Cartesian 
      * atom coordinates to the archive. In either case, it does not 
      * modify any Atom position values, or change the coordinates 
      * system setting of the AtomStorage.
      *
      * Usage: Atomic coordinates are Cartesian on input and output 
      * when this function is used by an Integrator to write a restart 
      * file, but are scaled / generalized on input and ouput when it 
      * is called by writeConfig() to implement the WRITE_CONFIG command 
      * within the Simulation::readCommand() function. 
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
