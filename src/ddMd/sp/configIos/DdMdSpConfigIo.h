#ifndef DDMD_SP_DDMD_CONFIG_IO_H
#define DDMD_SP_DDMD_CONFIG_IO_H

/*
* Simpatico - SpConfiguration Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/sp/configIos/SpConfigIo.h>
#include <ddMd/sp/storage/SpGroupStorage.h>
#include <ddMd/sp/storage/SpConfiguration.h>

#include <util/format/Int.h>

namespace DdMd
{

   class SpConfiguration;

   using namespace Util;

   /**
   * Native / default DdMd format for configuration files.
   *
   * \ingroup DdMd_SpConfigIo_Module
   */
   class DdMdSpConfigIo  : public SpConfigIo
   {

   public:

      /**
      * Default constructor.
      */
      DdMdSpConfigIo(bool hasMolecules = false);

      /**
      * Constructor.
      *
      * \param configuration parent SpConfiguration object.
      * \param hasMolecules true if file format has DdMd::AtomContext info
      */
      DdMdSpConfigIo(SpConfiguration& configuration, bool hasMolecules = false);

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
                     const char* nGroupLabel, SpGroupStorage<N>& groups);

      template <int N>
      int writeGroups(std::ofstream& file, const char* sectionLabel, 
                      const char* nGroupLabel, SpGroupStorage<N>& groups);

   };

   // Member functions templates

   /*
   * Private method to read SpGroup<N> objects.
   */
   template <int N>
   int DdMdSpConfigIo::readGroups(std::ifstream& file, 
                  const char* sectionLabel,
                  const char* nGroupLabel,
                  SpGroupStorage<N>& groups)
   {
      int nGroup;  // Total number of groups in file
      file >> Label(sectionLabel);
      file >> Label(nGroupLabel) >> nGroup;
      SpGroup<N>* groupPtr;
      for (int i = 0; i < nGroup; ++i) {
         groupPtr = groups.newPtr();
         file >> *groupPtr;
      }
      return nGroup;
   }

   /*
   * Private method to write SpGroup<N> objects.
   */
   template <int N>
   int DdMdSpConfigIo::writeGroups(std::ofstream& file, const char* sectionLabel,
                  const char* nGroupLabel, SpGroupStorage<N>& groups)
   {
      SpConfiguration::BondIterator iter;
      int nGroup = configuration().bonds().size();

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
