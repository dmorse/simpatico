#ifndef UTIL_XML_END_TAG_CPP
#define UTIL_XML_END_TAG_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "XmlEndTag.h"
#include "ioUtil.h"

namespace Util
{

   XmlEndTag::XmlEndTag()
   {}

   bool XmlEndTag::match(const std::string& line, int begin)
   {
      label_ = "";
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
      next();
      if (isEnd()) return false;
      if (c() != '/') {
         //std::cout << "Missing slash" << std::endl;
         return false;
      }

      // Read label
      int beginLabel, endLabel;
      next();
      skip();
      if (isEnd()) return false;
      beginLabel = cursor();
      while (c() != '>' && c() != ' ') {
         next();
         if (isEnd()) return false;
      }
      endLabel = cursor();
      label_ = string().substr(beginLabel, endLabel - beginLabel);
    
      skip();
      if (isEnd()) return false;
      if (c() == '>') {
         return true;
      } else {
         return false;
      }

   }

}
#endif
