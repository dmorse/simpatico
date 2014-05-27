#ifndef MDCF_DDMD_CONFIG_IO_H
#define MDCF_DDMD_CONFIG_IO_H

/*
* Simpatico - Storage Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mdCf/configIos/ConfigIo.h>
#include <mdCf/storage/GroupStorage.h>
#include <mdCf/storage/Storage.h>

#include <util/format/Int.h>

namespace MdCf
{

   class Storage;

   using namespace Util;

   /**
   * Native / default DdMd format for configuration files.
   *
   * \ingroup MdCf_ConfigIo_Module
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
      * \param storage parent Storage object.
      * \param hasMolecules true if file format has DdMd::AtomContext info
      */
      DdMdConfigIo(Storage& storage, bool hasMolecules = false);

      /**
      * Read configuration file in DdMd default format.
      *
      * \param file input file stream
      */
      virtual void readConfig(std::ifstream& file);

      /**
      * Write configuration file in DdMd default format.
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
      Storage::BondIterator iter;
      int nGroup = storage().bonds().size();

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
