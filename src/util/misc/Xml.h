#ifndef UTIL_XML_H
#define UTIL_XML_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ParserString.h"
#include <sstream>
#include <vector>

namespace Util
{





   /**
   * Parser for an XML attribute.
   */
   class XmlAttribute : public ParserString
   {

   public:

      /**
      * Constructor
      */ 
      XmlAttribute();

      /**
      * Return true if an attribute is found, false otherwise.
      */
      bool match(const std::string& string, int begin);

      /**
      * If successful return true and advance cursor or parent parser.
      *
      * \param parser parent parser object
      */
      bool match(ParserString& parser);

      /**
      * Return label string.
      */
      const std::string& label()
      {  return label_; }

      /**
      * Return value string, without quotes.
      */
      std::stringstream& value()
      {  return value_; }

   private:

      std::string label_;
      std::stringstream value_;

   };

   /**
   * Parser for an XML start tag.
   * 
   * Usage:
   * \code
   *    XmlStartTag tag;
   *    std::string line;
   *    tag.matchLabel(line, 0);
   *    while (matchAttribute()) {
   *       // process tag.attribute();
   *    } 
   *    if (!tag.isEnd()) {
   *       UTIL_THROW("Ill-formed start tag");
   *    } 
   * \endcode
   */
   class XmlStartTag : public ParserString
   {

   public:

      /**
      * Constructor
      */ 
      XmlStartTag();

      /**
      * Return true if an attribute is found, false otherwise.
      *
      * \param string containing text of XML tag
      * \param begin  index of first character
      */
      bool matchLabel(const std::string& string, int begin);

      bool matchAttribute();

      /**
      * Return attribute.
      */
      XmlAttribute& attribute()
      {  return attribute_; }

      bool isEnd()
      {  return isEnd_; }

   private:

      XmlAttribute attribute_;

      std::string label_;

      bool isEnd_;

   };

}
#endif
