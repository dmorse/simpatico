#ifndef SPAN_HOOMD_CONFIG_READER_H
#define SPAN_HOOMD_CONFIG_READER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <spAn/config/ConfigReader.h>   // base class
#include <spAn/config/TypeMap.h>        // member

#include <iostream>

namespace Util {
   class XmlStartTag;
}

namespace SpAn
{

   class Configuration;
   template <int N> class GroupStorage;

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

      /**
      * Read auxiliary file with type map information.
      *
      * \param file input file stream
      */
      virtual void readAuxiliaryFile(std::ifstream& file);

   private:

      TypeMap atomTypeMap_;
      TypeMap bondTypeMap_;
      TypeMap angleTypeMap_;
      TypeMap dihedralTypeMap_;

      /**
      * Read box data node.
      *
      * \param start XmlStartTag, after node name has been read
      * \param file  configuration file
      */
      void readBox(Util::XmlStartTag& start, std::istream& file);

      /**
      * Proces "num" attribute of atom or group node Xml start tag.
      *
      * \param nAtom number of atoms currently in AtomStorage
      * \return number of items expected in data node
      */
      int readNumberAttribute(Util::XmlStartTag& start, int nAtom);

      /**
      * Check xml end tag, throw exception if no match.
      * 
      * \param file input file
      * \param name expected name of end tag
      */
      void endTag(std::istream& file, const std::string& name);

      /**
      * Read position data node.
      *
      * \param start XmlStartTag, after node name has been read
      * \param file  configuration file
      */
      void readPosition(Util::XmlStartTag& start, std::istream& file);

      /**
      * Read velocity data node.
      *
      * \param start XmlStartTag, after node name has been read
      * \param file  configuration file
      */
      void readVelocity(Util::XmlStartTag& start, std::istream& file);

      /**
      * Read atom type data node.
      *
      * \param start XmlStartTag, after node name has been read
      * \param file  configuration file
      */
      void readType(Util::XmlStartTag& start, std::istream& file);

      /**
      * Read and discard and atom type data node.
      *
      * \param start XmlStartTag, after node name has been read
      * \param file  configuration file
      */
      void readAtomIgnore(Util::XmlStartTag& start, std::istream& file);

      /**
      * Read bond data node.
      *
      * \param start XmlStartTag, after node name has been read
      * \param file  configuration file
      */
      void readBond(Util::XmlStartTag& start, std::istream& file);

      /**
      * Read bond data node.
      *
      * \param start XmlStartTag, after node name has been read
      * \param file configuration file
      */
      template <int N>
      void readGroups(Util::XmlStartTag& start, std::istream& file,
                      GroupStorage<N>& storage, const TypeMap& map);

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
