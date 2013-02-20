#ifndef DDMD_DDMD_ORDERED_CONFIG_IO_H
#define DDMD_DDMD_ORDERED_CONFIG_IO_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/configIos/ConfigIo.h>
#include <util/space/Vector.h>

#include <vector>

namespace DdMd
{

   class Simulation;

   using namespace Util;

   /**
   * Native / default format for configuration files.
   *
   * DdMdOrderedConfigIo is a ConfigIo that implements the default
   * configuration file format for ddSim.
   *
   * \ingroup DdMd_ConfigIo_Module
   */
   class DdMdOrderedConfigIo  : public ConfigIo
   {

   public:

      /**
      * Default constructor.
      */
      DdMdOrderedConfigIo();

      /**
      * Constructor.
      *
      * \param simulation parent Simulation object.
      */
      DdMdOrderedConfigIo(Simulation& simulation);

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
      int readGroups(std::istream& file, 
                     const char* sectionLabel, const char* nGroupLabel,
                     GroupDistributor<N>& distributor);

      /**
      * Write Group<N> objects to file. 
      */
      template <int N>
      int writeGroups(std::ostream& file, 
                      const char* sectionLabel, const char* nGroupLabel,
                      GroupStorage<N>& storage, GroupCollector<N>& collector);
   
   };

}
#endif
