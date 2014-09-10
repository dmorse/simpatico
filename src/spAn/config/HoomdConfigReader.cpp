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
      if (!start.matchLabel(line, 0)) {
         UTIL_THROW("Missing hoomd_xml start tag");
      }
      if (start.label() != "hoomd_xml") {
         Log::file() << line << std::endl;
         Log::file() << "Start label = " << start.label() << std::endl;
         UTIL_THROW("Incorrect hoomd_xml start tag label");
      }
      XmlAttribute attribute;
      while (start.matchAttribute(attribute)) {
         UTIL_THROW("Found attribute of hoomd_xml tag");
      }
      if (!start.endBracket()) {
          Log::file() << line << std::endl;
          UTIL_THROW("Missing end bracket for hoomd_xml tag");
      }

      // Read data nodes
      XmlEndTag end;
      std::string name;
      bool finish = false;
      while (!finish) {
         
         getNextLine(file, line);

         // Attempt to match hoomd_xml end tag
         if (end.match(line, 0)) {

            if (end.label() != "hoomd_xml") {
               Log::file() << "End label = " << end.label() << std::endl;
               UTIL_THROW("Incorrect hoomd_xml end tag label");
            }
            finish = true;

         } else // Attempt to match a data node
         if (start.matchLabel(line, 0)) {

            name = start.label();
            Log::file() << std::endl;
            Log::file() << "node name = " << start.label() << std::endl;

            // Select function to read the remainder of the node, including attributes
            if (name == "atom") {
               readAtomNode(start, file);
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

   }
   
   void HoomdConfigReader::readAtomNode(Util::XmlStartTag& start, std::istream& file)
   {
       // Process data node start tag
       XmlAttribute attribute;
       while (start.matchAttribute(attribute)) {
           Log::file() << attribute.label() << std::endl;
       }
       if (!start.endBracket()) {
          UTIL_THROW("Missing end bracket");
       }

       // Process body (e.g., atom positions

       // Process data node end tag
       std::string line;
       XmlEndTag end;
       getNextLine(file, line);
       if (!end.match(line, 0)) {
          UTIL_THROW("Missing end tag for data node");
       }
       if (end.label() != start.label()) {
          UTIL_THROW("End tag label != start tag label");
       }
        
   }

      #if 0
      // Read and broadcast boundary
      file >> Label("BOUNDARY");
      file >> configuration().boundary();

      // Read and distribute atoms

      // Read atoms
      Atom* atomPtr;
      // atomCapacity = maximum allowed id + 1
      int atomCapacity = configuration().atoms().capacity(); 
      int nAtom;          
      file >> Label("ATOMS");
      file >> Label("nAtom") >> nAtom;
      for (int i = 0; i < nAtom; ++i) {

         // Get pointer to new atom 
         atomPtr = configuration().atoms().newPtr();
 
         file >> atomPtr->id;
         if (atomPtr->id < 0) {
            std::cout << "atom id =" << atomPtr->id << std::endl;
            UTIL_THROW("Negative atom id");
         }
         if (atomPtr->id >= atomCapacity) {
            std::cout << "atom id      =" << atomPtr->id << std::endl;
            std::cout << "atomCapacity =" << atomCapacity << std::endl;
            UTIL_THROW("Invalid atom id");
         }
         file >> atomPtr->typeId;
         if (hasMolecules_) {
            file >> atomPtr->speciesId;
            if (atomPtr->speciesId < 0) {
               std::cout << "species Id  =" << atomPtr->speciesId << std::endl;
               UTIL_THROW("Negative species id");
            }
            file >> atomPtr->moleculeId; 
            if (atomPtr->moleculeId < 0) {
               std::cout << "molecule Id =" << atomPtr->moleculeId << std::endl;
               UTIL_THROW("Negative molecule id");
            }
            file >> atomPtr->atomId;
            if (atomPtr->atomId < 0) {
               std::cout << "atom id     =" << atomPtr->atomId << std::endl;
               UTIL_THROW("Negative atom id in molecule");
            }
         }
         file >> atomPtr->position;
         file >> atomPtr->velocity;

         // Finalize addition of new atom
         configuration().atoms().add();
      }

      // Read Covalent Groups
      #ifdef INTER_BOND
      if (configuration().bonds().capacity()) {
         readGroups(file, "BONDS", "nBond", configuration().bonds());
         //if (maskPolicy == MaskBonded) {
         //   setAtomMasks();
         //}
      }
      #endif

      #ifdef INTER_ANGLE
      if (configuration().angles().capacity()) {
         readGroups(file, "ANGLES", "nAngle", configuration().angles());
      }
      #endif

      #ifdef INTER_DIHEDRAL
      if (configuration().dihedrals().capacity()) {
         readGroups(file, "DIHEDRALS", "nDihedral", configuration().dihedrals());
      }
      #endif

      #endif // if 0

}
#endif
