#ifndef SPAN_HOOMD_CONFIG_READER_CPP
#define SPAN_HOOMD_CONFIG_READER_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "HoomdConfigReader.h"

#include <spAn/chemistry/Atom.h>
#include <spAn/chemistry/Group.h>
#include <spAn/chemistry/Species.h>
//#include <spAn/chemistry/MaskPolicy.h>
#include <spAn/storage/Configuration.h>

#include <util/space/Vector.h>
#include <util/xmltag/XmlXmlTag.h>
#include <util/xmltag/XmlStartTag.h>
#include <util/xmltag/XmlAttribute.h>
#include <util/xmltag/XmlEndTag.h>
#include <util/misc/ioUtil.h>

namespace SpAn
{

   using namespace Util;

   /*
   * Constructor.
   */
   HoomdConfigReader::HoomdConfigReader(Configuration& configuration)
    : ConfigReader(configuration)
   {  setClassName("HoomdConfigReader"); }

   /*
   * Read a configuration file.
   */
   void HoomdConfigReader::readConfig(std::ifstream& file)
   {
      // Precondition
      if (!file.is_open()) {  
         UTIL_THROW("Error: File is not open"); 
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
      int dimension;
      double vizsigma;
      while (start.matchAttribute(attribute)) {
         if (attribute.label() == "time_step") {
            attribute.value() >> timestep;
         } else 
         if (attribute.label() == "dimension") {
            attribute.value() >> dimension;
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
            Log::file() << std::endl;
            Log::file() << "node name = " << start.label() << std::endl;

            // Select node type
            if (name == "position") {
               readPosition(start, file);
            } else
            if (name == "velocity") {
               readVelocity(start, file);
            } else
            if (name == "box") {
               Log::file() << "Box node not yet implemented" << std::endl;
            } else {
               UTIL_THROW("Unknown node name");
            }

         } else {

            Log::file() << line << std::endl;
            UTIL_THROW("Error unrecognized node");

         }
         
      }

      getNextLine(file, line);
      end.match("hoomd_xml", line, 0); 
   }





   void HoomdConfigReader::readBox(Util::XmlStartTag& start, 
                                        std::istream& file)
   {
      // Process data node start tag
      XmlAttribute attribute;
      Vector lengths;
      while (start.matchAttribute(attribute)) {
         if (attribute.label() == "lx") {
            attribute.value() >> lengths[0];
         } else 
         if (attribute.label() == "ly") {
            attribute.value() >> lengths[1];
         } else 
         if (attribute.label() == "lz") {
            attribute.value() >> lengths[2];
         } else {
            Log::file() << attribute.label() << std::endl;
            UTIL_THROW("Unknown attribute");
         }
      }
      start.finish();

      configuration().boundary().setOrthorhombic(lengths);

   }

   void HoomdConfigReader::readPosition(Util::XmlStartTag& start, 
                                        std::istream& file)
   {
      // Count atoms in storage
      AtomStorage& storage = configuration().atoms();
      int n = storage.size();
      bool hasAtoms;
      if (n == 0) {
         hasAtoms = false; 
      } else {
         hasAtoms = true; 
      }

      // Process data node start tag
      XmlAttribute attribute;
      while (start.matchAttribute(attribute)) {
         if (attribute.label() == "num") {
            int num;
            attribute.value() >> num;
            if (hasAtoms && num != n) {
               UTIL_THROW("Inconsistent number of atoms");
            } 
         } else {
            Log::file() << attribute.label() << std::endl;
            UTIL_THROW("Unknown attribute");
         }
      }
      start.finish();

      // Read atom positions
      Atom* atomPtr;
      if (hasAtoms) {
         for (int i = 0; i < n; ++i) {
            atomPtr = storage.ptr(i);
            if (!atomPtr) {
               UTIL_THROW("Atom not found");   
            }
            file >> atomPtr->position;
         }
      } else {
         for (int i = 0; i < n; ++i) {
            atomPtr = storage.newPtr();
            atomPtr->id = i;
            file >> atomPtr->position;
            storage.add();
         }
      }

      // Process end tag
      std::string line;
      XmlEndTag end;
      getNextLine(file, line);
      end.match("position", line, 0);
   }

   void HoomdConfigReader::readVelocity(Util::XmlStartTag& start, 
                                        std::istream& file)
   {
      // Count atoms in storage
      AtomStorage& storage = configuration().atoms();
      int n = storage.size();
      bool hasAtoms;
      if (n == 0) {
         hasAtoms = false; 
      } else {
         hasAtoms = true; 
      }

      // Process data node start tag
      XmlAttribute attribute;
      while (start.matchAttribute(attribute)) {
         if (attribute.label() == "num") {
            int num;
            attribute.value() >> num;
            if (hasAtoms && num != n) {
               UTIL_THROW("Inconsistent number of atoms");
            } 
         } else {
            Log::file() << attribute.label() << std::endl;
            UTIL_THROW("Unknown attribute");
         }
      }
      start.finish();

      // Read atom velocitys
      Atom* atomPtr;
      if (hasAtoms) {
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

      // Process end tag
      std::string line;
      XmlEndTag end;
      getNextLine(file, line);
      end.match("velocity", line, 0);
   }

   void HoomdConfigReader::readType(Util::XmlStartTag& start, 
                                    std::istream& file)
   {
      // Count atoms in storage
      AtomStorage& storage = configuration().atoms();
      int n = storage.size();
      bool hasAtoms;
      if (n == 0) {
         hasAtoms = false; 
      } else {
         hasAtoms = true;
      }

      // Process data node start tag
      XmlAttribute attribute;
      while (start.matchAttribute(attribute)) {
         if (attribute.label() == "num") {
            int num;
            attribute.value() >> num;
            if (hasAtoms && num != n) {
               UTIL_THROW("Inconsistent number of atoms");
            } 
         } else {
            Log::file() << attribute.label() << std::endl;
            UTIL_THROW("Unknown attribute");
         }
      }
      start.finish();

      // Read atom velocities
      Atom* atomPtr;
      if (hasAtoms) {
         for (int i = 0; i < n; ++i) {
            atomPtr = storage.ptr(i);
            if (!atomPtr) {
               UTIL_THROW("Atom not found");   
            }
            file >> atomPtr->typeId;
         }
      } else {
         for (int i = 0; i < n; ++i) {
            atomPtr = storage.newPtr();
            atomPtr->id = i;
            file >> atomPtr->typeId;
            storage.add();
         }
      }

      // Process end tag
      std::string line;
      XmlEndTag end;
      getNextLine(file, line);
      end.match("type", line, 0);
   }

}
#endif
