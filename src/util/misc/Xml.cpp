#ifndef UTIL_XML_CPP
#define UTIL_XML_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Xml.h"
#include "ioUtil.h"

namespace Util
{

   XmlAttribute::XmlAttribute()
    : 
      ParserString(),
      label_(""),
      value_()
   {
      value_.str("");
      value_.clear(); 
   }

   bool XmlAttribute::match(const std::string& line, int begin)
   {
      setString(line, begin);
      label_ = "";
      value_.str("");
      value_.clear();

      // Read label
      int beginLabel, endLabel;
      skip();
      if (isEnd()) return false;
      if (c() == '=') return false;
      if (c() == '>') return false;
      beginLabel = cursor();
      while (c() != '=') {
         next();
         if (isEnd()) return false;
         if (c() == '>') return false;
      }
      endLabel = cursor();

      // Advance and skip white space
      next();
      skip();
      if (isEnd()) return false;

      // Read value
      int beginValue, endValue;
      char quote;
      if (c() == '\'' || c() == '\"') {
         quote = c();
      }
      next();
      if (isEnd()) return false;
      beginValue = cursor();
      while (c() != quote) {
         next();
         if (isEnd()) return false;
      }
      endValue = cursor();
      next();
  
      label_ = string().substr(beginLabel, endLabel - beginLabel); 
      rStrip(label_);
      value_.str(string().substr(beginValue, endValue - beginValue)); 
 
      return true;
   }

   bool XmlAttribute::match(ParserString& parser)
   {
      bool status = match(parser.string(), parser.cursor());
      if (status) {
         parser.setCursor(cursor());
      }
      return status;
   }

   XmlStartTag::XmlStartTag()
    : isEnd_(false)
   {}

   bool XmlStartTag::matchLabel(const std::string& line, int begin)
   {
      setString(line, begin);
      skip();
      if (isEnd()) return false;

      // Read opening bracket
      if (c() != '<') return false;

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
      skip();
      if (isEnd()) return false;

   }
         
   bool XmlStartTag::matchAttribute()
   {
      if (c() == '>') {
         isEnd_ = true;
         return false;
      } else {
         return attribute_.match(*this);
      }
   }

}
#endif
