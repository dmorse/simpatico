/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "HoomdConfigReader.h"

#include <tools/chemistry/Atom.h>
#include <tools/chemistry/Group.h>
#include <tools/chemistry/Species.h>
//#include <tools/chemistry/MaskPolicy.h>
#include <tools/storage/Configuration.h>

#include <util/space/Vector.h>
#include <util/xmltag/XmlXmlTag.h>
#include <util/xmltag/XmlStartTag.h>
#include <util/xmltag/XmlAttribute.h>
#include <util/xmltag/XmlEndTag.h>
#include <util/misc/ioUtil.h>

namespace Tools
{

   using namespace Util;

   /*
   * Constructor.
   */
   HoomdConfigReader::HoomdConfigReader(Configuration& configuration)
    : ConfigReader(configuration, true),
      hasTypeMaps_(false)
   {  setClassName("HoomdConfigReader"); }

   /*
   * Read auxiliary type map file.
   */
   void HoomdConfigReader::readAuxiliaryFile(std::ifstream& file)
   {
      bool notEnd;
      std::stringstream line;

      notEnd = getNextLine(file, line);
      if (notEnd) {
         checkString(line, "ATOM");
         checkString(line, "TYPES:");
         atomTypeMap_.read(file);
      }

      notEnd = getNextLine(file, line);
      if (notEnd) {
         checkString(line, "BOND");
         checkString(line, "TYPES:");
         bondTypeMap_.read(file);
      }

      notEnd = getNextLine(file, line);
      if (notEnd) {
         checkString(line, "ANGLE");
         checkString(line, "TYPES:");
         angleTypeMap_.read(file);
      }

      hasTypeMaps_ = true;
   }

   /*
   * Read a configuration file.
   */
   void HoomdConfigReader::readConfig(std::ifstream& file)
   {
      // Preconditions
      if (!file.is_open()) {  
         UTIL_THROW("Error: File is not open"); 
      }
      if (!hasTypeMaps_) {
         UTIL_THROW("Error: A type map file must be read before config"); 
      }

      // Set flags
     
      std::string line;

      // Read XML tag
      getNextLine(file, line);
      XmlXmlTag xmlTag;
      if (!xmlTag.match(line, 0)) {
         UTIL_THROW("Missing XML tag");
      }

      // Read hoom_xml start tag
      getNextLine(file, line);
      XmlStartTag start;
      start.matchLabel("hoomd_xml", line, 0); 
      XmlAttribute attribute;
      std::string version;
      while (start.matchAttribute(attribute)) {
         if (attribute.label() == "version") {
            attribute.value() >> version;
         } else {
           Log::file() << "attribute = " << attribute.label() << std::endl;
           UTIL_THROW("Unknown attribute of hoomd_xml tag");
         }
      }

      // Read configuration start tag
      getNextLine(file, line);
      start.matchLabel("configuration", line, 0); 
      int timestep;
      int dimensions;
      int natoms;
      double vizsigma;
      while (start.matchAttribute(attribute)) {
         if (attribute.label() == "time_step") {
            attribute.value() >> timestep;
         } else 
         if (attribute.label() == "dimensions") {
            attribute.value() >> dimensions;
            if (dimensions != 3) {
               UTIL_THROW("hoomd_xml dimensions attribute not equal to 3");
            }
         } else
         if (attribute.label() == "natoms") {
            attribute.value() >> natoms;
         } else 
         if (attribute.label() == "vizsigma") {
            attribute.value() >> vizsigma;
         } else {
            Log::file() << "attribute = " << attribute.label() << std::endl;
            UTIL_THROW("Unknown attribute of hoomd_xml tag");
         }
      }
      start.finish();

      // Read data nodes
      XmlEndTag end;
      std::string name;
      bool finish = false;
      while (!finish) {
         
         getNextLine(file, line);

         // Attempt to match hoomd_xml end tag
         if (end.match(line, 0)) {

            if (end.label() != "configuration") {
               Log::file() << "End label = " << end.label() << std::endl;
               UTIL_THROW("Incorrect configuration end tag label");
            }
            finish = true;

         } else // Attempt to match a data node
         if (start.matchLabel(line, 0)) {

            name = start.label();
            // Log::file() << std::endl;
            // Log::file() << "node name = " << start.label() << std::endl;

            // Select node type
            if (name == "box") {
               readBox(start, file);
            } else
            if (name == "position") {
               readPosition(start, file);
            } else
            if (name == "velocity") {
               readVelocity(start, file);
            } else 
            if (name == "type") {
               readType(start, file);
            } else 
            if (name == "mass") {
               readAtomIgnore(start, file);
            } else 
            if (name == "charge") {
               readAtomIgnore(start, file);
            } else 
            if (name == "diameter") {
               readAtomIgnore(start, file);
            } else 
            if (name == "body") {
               readAtomIgnore(start, file);
            } else 
            if (name == "bond") {
               readGroups<2>(start, file, configuration().bonds(), 
                             bondTypeMap_);
            }
            #ifdef SIMP_ANGLE
            else if (name == "angle") {
               readGroups<3>(start, file, configuration().angles(), 
                             angleTypeMap_);
            } 
            #endif
            #ifdef SIMP_DIHEDRAL
            else if (name == "dihedral") {
               readGroups<4>(start, file, configuration().dihedrals(), 
                             dihedralTypeMap_);
            } else if (name == "improper") {
               readGroups<4>(start, file, configuration().impropers(), 
                             improperTypeMap_);
            } 
            #endif
            else {
               UTIL_THROW("Unknown node name");
            }

         } else {

            Log::file() << line << std::endl;
            UTIL_THROW("Expected data node start tag or configuration end tag");

         }
         
      }

      // Read hoomd_xml end line
      getNextLine(file, line);
      end.match("hoomd_xml", line, 0); 

      // If species are declared, attempt to set atom context info,
      // assuming that atom ids are ordered by molecule and species.

      if (configuration().nSpecies() > 0) {
         bool success;
         success = setAtomContexts();
         if (success) {
            addAtomsToSpecies();
         }
      }

   }





   void HoomdConfigReader::readBox(Util::XmlStartTag& start, 
                                        std::istream& file)
   {
      // Process data node start tag
      XmlAttribute attribute;
      Vector lengths;
      double xy, xz, yz;
      while (start.matchAttribute(attribute)) {
         if (attribute.label() == "lx") {
            attribute.value() >> lengths[0];
         } else 
         if (attribute.label() == "ly") {
            attribute.value() >> lengths[1];
         } else 
         if (attribute.label() == "lz") {
            attribute.value() >> lengths[2];
         } else 
         if (attribute.label() == "xy") {
            attribute.value() >> xy;
         } else 
         if (attribute.label() == "xz") {
            attribute.value() >> xz;
         } else 
         if (attribute.label() == "yz") {
            attribute.value() >> yz;
         } else {
            Log::file() << attribute.label() << std::endl;
            UTIL_THROW("Unknown attribute");
         }
      }
      start.finish();

      configuration().boundary().setOrthorhombic(lengths);

   }

   /*
   * Finish processing start tag.
   */
   int HoomdConfigReader::readNumberAttribute(Util::XmlStartTag& start, int nAtom)
   {
      XmlAttribute attribute;
      while (start.matchAttribute(attribute)) {
         if (attribute.label() == "num") {
            int num;
            attribute.value() >> num;
            if (nAtom > 0 && num != nAtom) {
               UTIL_THROW("Inconsistent number of atoms");
            } 
            if (nAtom == 0) {
               nAtom = num;
            }
         } else {
            Log::file() << attribute.label() << std::endl;
            UTIL_THROW("Unknown attribute");
         }
      }
      start.finish();
      return nAtom;
   }

   /*
   * Check end tag.
   */
   void HoomdConfigReader::endTag(std::istream& file, const std::string& name)
   {
      std::string line;
      XmlEndTag end;
      getNextLine(file, line);
      end.match(name, line, 0);
   }

   /*
   * Read position data node.
   */
   void HoomdConfigReader::readPosition(Util::XmlStartTag& start, 
                                        std::istream& file)
   {
      AtomStorage& storage = configuration().atoms();
      Boundary& boundary = configuration().boundary();
      int nAtom = storage.size();
      int n = readNumberAttribute(start, nAtom);

      Atom* atomPtr;
      if (nAtom > 0) {
         for (int i = 0; i < n; ++i) {
            atomPtr = storage.ptr(i);
            if (!atomPtr) {
               UTIL_THROW("Atom not found");   
            }
            file >> atomPtr->position;
            boundary.shift(atomPtr->position);
         }
      } else {
         for (int i = 0; i < n; ++i) {
            atomPtr = storage.newPtr();
            atomPtr->id = i;
            file >> atomPtr->position;
            boundary.shift(atomPtr->position);
            storage.add();
         }
      }
      endTag(file, "position");
   }

   /*
   * Read velocity data node.
   */
   void HoomdConfigReader::readVelocity(Util::XmlStartTag& start, 
                                        std::istream& file)
   {
      AtomStorage& storage = configuration().atoms();
      int nAtom = storage.size();
      int n = readNumberAttribute(start, nAtom);

      Atom* atomPtr;
      if (nAtom > 0) {
         for (int i = 0; i < n; ++i) {
            atomPtr = storage.ptr(i);
            if (!atomPtr) {
               UTIL_THROW("Atom not found");   
            }
            file >> atomPtr->velocity;
         }
      } else {
         for (int i = 0; i < n; ++i) {
            atomPtr = storage.newPtr();
            atomPtr->id = i;
            file >> atomPtr->velocity;
            storage.add();
         }
      }
      endTag(file, "velocity");
   }

   /*
   * Read type data node.
   */
   void HoomdConfigReader::readType(Util::XmlStartTag& start, 
                                    std::istream& file)
   {
      AtomStorage& storage = configuration().atoms();
      int nAtom = storage.size();
      int n = readNumberAttribute(start, nAtom);
      int typeId;
      std::string typeName;

      Atom* atomPtr;
      if (nAtom > 0) {
         for (int i = 0; i < n; ++i) {
            atomPtr = storage.ptr(i);
            if (!atomPtr) {
               UTIL_THROW("Atom not found");   
            }
            file >> typeName;
            typeId = atomTypeMap_.id(typeName);
            atomPtr->typeId = typeId;
         }
      } else {
         for (int i = 0; i < n; ++i) {
            atomPtr = storage.newPtr();
            atomPtr->id = i;
            file >> typeName;
            typeId = atomTypeMap_.id(typeName);
            atomPtr->typeId = typeId;
            storage.add();
         }
      }
      endTag(file, "type");
   }

   /*
   * Read and discard atom data 
   */
   void HoomdConfigReader::readAtomIgnore(Util::XmlStartTag& start, 
                                          std::istream& file)
   {
      AtomStorage& storage = configuration().atoms();
      int nAtom = storage.size();
      int n = readNumberAttribute(start, nAtom);

      std::string dummy;
      for (int i = 0; i < n; ++i) {
         file >> dummy;
      }

      endTag(file, start.label());
      #if 0
      XmlEndTag end;
      std::string line;
      getNextLine(file, line);
      if (!end.match(line, 0)) {
         UTIL_THROW("No end tag");
      }
      #endif
   }

   /*
   * Read type data node.
   */
   void HoomdConfigReader::readBond(Util::XmlStartTag& start, 
                                    std::istream& file)
   {
      GroupStorage < 2> & storage = configuration().bonds();
      int nGroup = storage.size();
      if (nGroup > 0) {
         UTIL_THROW("Group storage not empty");
      }
      int n = readNumberAttribute(start, nGroup);
      int typeId;
      std::string typeName;

      int i, j;
      Group<2>* groupPtr;
      for (i = 0; i < n; ++i) {
         groupPtr = storage.newPtr();
         groupPtr->id = i;
         file >> typeName;
         typeId = bondTypeMap_.id(typeName);
         groupPtr->typeId = typeId;
         for (j = 0; j < 2; ++j) {
            file >> groupPtr->atomIds[j];
         }
      }
      endTag(file, start.label());
   }

   /*
   * Read group data node.
   */
   template <int N>
   void HoomdConfigReader::readGroups(Util::XmlStartTag& start, 
                                      std::istream& file,
                                      GroupStorage<N>& storage,
                                      const TypeMap& map)
   {
      int nGroup = storage.size();
      if (nGroup > 0) {
         UTIL_THROW("Group storage not empty");
      }
      int n = readNumberAttribute(start, nGroup);
      int typeId;
      std::string typeName;

      if (n > 0) {
         int i, j;
         Group<N>* groupPtr;
         for (i = 0; i < n; ++i) {
            groupPtr = storage.newPtr();
            groupPtr->id = i;
            file >> typeName;
            typeId = map.id(typeName);
            groupPtr->typeId = typeId;
            for (j = 0; j < 2; ++j) {
               file >> groupPtr->atomIds[j];
            }
         }
      }
      endTag(file, start.label());
   }

}
