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
   * Native / default format for configuration files.
   *
   * LammpsConfigIo is a ConfigIo that implements the default
   * configuration file format for ddSim.
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
      * This routine opens and reads a file on the master, and distributes
      * atom data among the processors.
      *
      * \param file input file stream
      * \param maskPolicy MaskPolicy to be used in setting atom masks
      */
      virtual void readConfig(std::istream& file, MaskPolicy maskPolicy);

      /**
      * Write configuration file.
      *
      * This routine opens and writes a file on the master,
      * collecting atom data from all processors.
      *
      * \param file output file stream
      */
      virtual void writeConfig(std::ostream& file);
   
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
      * Array of atoms, ordered by global index.
      */
      //int N;
      //std::vector<IoGroup <N> > groups_;

      /**
      * Read Group<N> objects from file. 
      */
      template <int N>
      void readGroups(std::istream& file, const char* sectionLabel, 
                      int nGroup, GroupDistributor<N>& distributor);

      /**
      * Write Group<N> objects to file. 
      */
      template <int N>
      void writeGroups(std::ostream& file, const char* sectionLabel,
                       GroupStorage<N>& storage, GroupCollector<N>& collector);
   
   };

}
#endif