#ifndef TOOLS_HOOMD_CONFIG_READER_H
#define TOOLS_HOOMD_CONFIG_READER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <tools/config/ConfigReader.h>   // base class
#include <tools/config/TypeMap.h>        // member

#include <iostream>

namespace Util {
   class XmlStartTag;
}

namespace Tools
{

   class Configuration;
   template <int N> class GroupStorage;

   using namespace Util;

   /**
   * Reader for Hoomd-blue XML configuration files.
   *
   * This config file reader requires access to an 
   * auxiliary file containing the mapping between
   * atom type strings used in the Hoomd XML file
   * format and the integer atom type ids used 
   * internally by simpatico. The required syntax
   * for the READ_CONFIG command in the mdPp command
   * file is thus
   * \code
   *    READ_CONFIG auxiliaryFile configFile
   * \endcode
   * where auxiliaryFile is the name of the auxilary
   * type map file and configFile is the name of the
   * new Hoomd XML configuration file. 
   *
   * \sa \ref tools_config_HoomdTypeMap_page
   *
   * \ingroup Tools_ConfigReader_Module
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
      TypeMap improperTypeMap_;
      bool hasTypeMaps_;

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
      * Read Group<N> (bond, angle, ...) data node.
      *
      * \param start XmlStartTag, after node label has been matched
      * \param file configuration file
      */
      template <int N>
      void readGroups(Util::XmlStartTag& start, std::istream& file,
                      GroupStorage<N>& storage, const TypeMap& map);

   };

}
#endif
