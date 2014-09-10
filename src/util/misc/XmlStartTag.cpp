#ifndef UTIL_XML_START_TAG_CPP
#define UTIL_XML_START_TAG_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "XmlStartTag.h"
#include "ioUtil.h"

namespace Util
{

   XmlStartTag::XmlStartTag()
    : endBracket_(false)
   {}

   bool XmlStartTag::matchLabel(const std::string& line, int begin)
   {
      label_ = "";
      endBracket_ = false;
      setString(line, begin);

      // Skip leading white space
      if (isEnd()) return false;
      skip();
      if (isEnd()) return false;

      // Read opening bracket
      if (c() != '<') {
         //std::cout << "Missing opening bracket" << std::endl;
         return false;
      }

      // Read label
      int beginLabel, endLabel;
      next();
      skip();
      if (isEnd()) return false;
      beginLabel = cursor();
      while (c() != ' ') {
         next();
         if (isEnd()) return false;
      }
      endLabel = cursor();
      label_ = string().substr(beginLabel, endLabel - beginLabel);
      return true;

   }
         
   bool XmlStartTag::matchAttribute(XmlAttribute& attribute)
   {
      skip();
      if (isEnd()) return false;
      if (c() == '>') {
         endBracket_ = true;
         return false;
      } else if (c() == '\/') {
         next();
         if (isEnd()) return false;
         if (c() == '>') {
            endBracket_ = true;
         }
         return false;
      } else {
         return attribute.match(*this);
      }
   }

}
#endif
