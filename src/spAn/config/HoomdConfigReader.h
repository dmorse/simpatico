#ifndef SPAN_HOOMD_CONFIG_READER_H
#define SPAN_HOOMD_CONFIG_READER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <spAn/config/ConfigReader.h>
//#include <spAn/storage/GroupStorage.h>
//#include <spAn/storage/Configuration.h>

#include <iostream>

namespace Util {
   class XmlStartTag;
}

namespace SpAn
{

   class Configuration;

   using namespace Util;

   /**
   * Reader for Hoomd-blue XML configuration files.
   *
   * \ingroup SpAn_ConfigReader_Module
   */
   class HoomdConfigReader  : public ConfigReader
   {

   public:

      /**
      * Default constructor.
      */
      HoomdConfigReader();

      /**
      * Constructor.
      *
      * \param configuration parent Configuration object.
      */
      HoomdConfigReader(Configuration& configuration);

      /**
      * Read configuration file in DdMd default format.
      *
      * \param file input file stream
      */
      virtual void readConfig(std::ifstream& file);

   private:

      #if 0
      template <int N>
      int readGroups(std::ifstream& file, const char* sectionLabel, 
                     const char* nGroupLabel, GroupStorage<N>& groups);
      #endif

      void processNode(Util::XmlStartTag& start, std::istream& file);

   };

   // Member functions templates

   #if 0
   /*
   * Private method to read Group<N> objects.
   */
   template <int N>
   int HoomdConfigReader::readGroups(std::ifstream& file, 
                  const char* sectionLabel,
                  const char* nGroupLabel,
                  GroupStorage<N>& groups)
   {
      int nGroup;  // Total number of groups in file
      //file >> Label(sectionLabel);
      //file >> Label(nGroupLabel) >> nGroup;
      Group<N>* groupPtr;
      for (int i = 0; i < nGroup; ++i) {
         groupPtr = groups.newPtr();
         //file >> *groupPtr;
      }
      return nGroup;
   }
   #endif

}
#endif
