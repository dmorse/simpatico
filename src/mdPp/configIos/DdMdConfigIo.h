#ifndef MDPP_DDMD_CONFIG_IO_H
#define MDPP_DDMD_CONFIG_IO_H

/*
* Simpatico - Processor Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mdPp/configIos/ConfigIo.h>
#include <mdPp/processor/GroupStorage.h>
#include <mdPp/processor/Processor.h>

#include <util/format/Int.h>

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

      template <int N>
      int readGroups(std::ifstream& file, const char* sectionLabel, 
                     const char* nGroupLabel, GroupStorage<N>& groups);

      template <int N>
      int writeGroups(std::ofstream& file, const char* sectionLabel, 
                      const char* nGroupLabel, GroupStorage<N>& groups);

   };

   // Member functions templates

   /*
   * Private method to read Group<N> objects.
   */
   template <int N>
   int DdMdConfigIo::readGroups(std::ifstream& file, 
                  const char* sectionLabel,
                  const char* nGroupLabel,
                  GroupStorage<N>& groups)
   {
      int nGroup;  // Total number of groups in file
      file >> Label(sectionLabel);
      file >> Label(nGroupLabel) >> nGroup;
      Group<N>* groupPtr;
      for (int i = 0; i < nGroup; ++i) {
         groupPtr = groups.newPtr();
         file >> *groupPtr;
      }
      return nGroup;
   }

   /*
   * Private method to write Group<N> objects.
   */
   template <int N>
   int DdMdConfigIo::writeGroups(std::ofstream& file, const char* sectionLabel,
                  const char* nGroupLabel, GroupStorage<N>& groups)
   {
      Processor::BondIterator iter;
      int nGroup = processor().bonds().size();

      file << std::endl;
      file << sectionLabel << std::endl;
      file << nGroupLabel << Int(nGroup, 10) << std::endl;
      for (groups.begin(iter); iter.notEnd(); ++iter) {
         file << *iter << std::endl;
      }
      return nGroup;
   }

}
#endif
