#ifndef DDMD_LAMMPS_CONFIG_IO_H
#define DDMD_LAMMPS_CONFIG_IO_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/configIos/ConfigIo.h>
#include <ddMd/chemistry/Group.h>
#include <util/space/Vector.h>

#include <vector>

namespace DdMd
{

   class Simulation;

   using namespace Util;

   /**
   * Lammps data file format for configuration files.
   *
   * LammpsConfigIo is a ConfigIo that reads and writes the LAMMPS data file
   * format for configuration files. Atoms are output with ids in sequential
   * order.
   *
   * \ingroup DdMd_ConfigIo_Module
   */
   class LammpsConfigIo  : public ConfigIo
   {

   public:

      /**
      * Default constructor.
      */
      LammpsConfigIo();

      /**
      * Constructor.
      *
      * \param simulation parent Simulation object.
      */
      LammpsConfigIo(Simulation& simulation);

      /**
      * Read configuration file.
      *
      * Call on all processors. File is used and must be open only on master.
      * This routine opens and reads a file on the master, and distributes
      * atom data among the processors.
      *
      * \pre  There are no atoms, ghosts, or groups.
      * \pre  AtomStorage is set for scaled / generalized coordinates
      *
      * \post atomic coordinates are scaled / generalized
      * \post there are no ghosts
      *
      * \param file input file stream
      * \param maskPolicy MaskPolicy to be used in setting atom masks
      */
      virtual void readConfig(std::ifstream& file, MaskPolicy maskPolicy);

      /**
      * Write configuration file.
      *
      * Call on all processors. File is used and must be open only on master.
      * This routine writes a file on the master, collecting atom and group
      * data from all processors.
      *
      * \param file output file stream (must be open on master)
      */
      virtual void writeConfig(std::ofstream& file);
   
   private:

      int nAtomType_;

      int nBondType_;

      int nAngleType_;

      int nDihedralType_;

      int nImproperType_;

      struct IoAtom {
         Vector position;
         Vector velocity;
         int    typeId;
         int    id;
      };

      /**
      * Array of atoms, ordered by global index.
      */
      std::vector<IoAtom> atoms_;

      template <int N>
      struct IoGroup {
         int    id;
         Group<N> group;
      };
  
      /**
      * Read Group<N> objects from file. 
      */
      template <int N>
      void readGroups(std::ifstream& file, const char* sectionLabel, 
                      int nGroup, GroupDistributor<N>& distributor);

      /**
      * Write Group<N> objects to file. 
      */
      template <int N>
      void writeGroups(std::ofstream& file, const char* sectionLabel,
                       GroupStorage<N>& storage, GroupCollector<N>& collector);
   
   };

}
#endif
