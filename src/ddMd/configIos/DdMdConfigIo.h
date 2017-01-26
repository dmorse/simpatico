#ifndef DDMD_DDMD_CONFIG_IO_H
#define DDMD_DDMD_CONFIG_IO_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
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
      *
      * \param hasMolecules desired value of hasMolecules flag
      */
      DdMdConfigIo(bool hasMolecules = false);

      /**
      * Constructor.
      *
      * \param simulation parent Simulation object
      * \param hasMolecules desired value of hasMolecules flag
      */
      DdMdConfigIo(Simulation& simulation, bool hasMolecules = false);

      /**
      * Read configuration file.
      *
      * This routine opens and reads a file on the master, and distributes
      * atom data among the processors.
      *
      * \pre  There are no atoms, ghosts, or groups.
      * \pre  AtomStorage is set for scaled / generalized coordinates
      *
      * \post Atomic coordinates are scaled / generalized
      * \post There are no ghosts
      *
      * \param file input file stream (must be open on master)
      * \param maskPolicy MaskPolicy to be used in setting atom masks
      */
      virtual void readConfig(std::ifstream& file, MaskPolicy maskPolicy);

      /**
      * Write configuration file.
      *
      * This routine opens and writes a file on the master, collecting 
      * atom data from all processors. Atomic positions are written in
      * Cartesian coordinates, independent of the coordinate system used
      * by the AtomStorage.
      *
      * \param file output file stream (must be open on master)
      */
      virtual void writeConfig(std::ofstream& file);
   
   private:

      /**
      * Include AtomContext info in file format, if available? 
      */
      bool hasMolecules_;

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
