#ifndef UTIL_XML_END_TAG_H
#define UTIL_XML_END_TAG_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "XmlBase.h"
#include "XmlAttribute.h"
#include <sstream>
#include <vector>

namespace Util
{
   
   /**
   * Parser for an XML end tag.
   *
   * \ingroup Misc_Module
   */
   class XmlEndTag : public XmlBase
   {

   public:

      /**
      * Constructor
      */ 
      XmlEndTag();

      /**
      * Return true if an attribute is found, false otherwise.
      *
      * \param string containing text of XML tag
      * \param begin  index of first character
      */
      bool match(const std::string& string, int begin);

      /**
      * Label string.
      */
      const std::string label()
      {  return label_; }

   private:

      /**
      * Label string (name of XML element).
      */
      std::string label_;

   };

}
#endif
