#ifndef UTIL_XML_START_TAG_H
#define UTIL_XML_START_TAG_H

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
   * Parser for an XML start tag.
   * 
   * Usage:
   * \code
   *    XmlStartTag tag;
   *    XmlAttribute attribute;
   *    std::string line;
   *    tag.matchLabel(line, 0);
   *    while (matchAttribute(attribute)) {
   *       // process attribute;
   *    } 
   *    if (!tag.endBracket()) {
   *       UTIL_THROW("No end-bracket - ill formed start tag");
   *    } 
   * \endcode
   *
   * \ingroup XmlTag_Module
   */
   class XmlStartTag : public XmlBase
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

      /**
      * Attempt to match an attribute.
      *
      * \param  attribute on return, matched attribute, if any
      * \return true if an attribute is found, false otherwise
      */
      bool matchAttribute(XmlAttribute& attribute);

      /**
      * Label string.
      */
      const std::string label()
      {  return label_; }

      /**
      * True if a closing bracket was found.
      */
      bool endBracket()
      {  return endBracket_; }

   private:

      /**
      * Label string (name of XML element).
      */
      std::string label_;

      /**
      * Set true when end bracket found, false until then.
      */
      bool endBracket_;

   };

}
#endif
