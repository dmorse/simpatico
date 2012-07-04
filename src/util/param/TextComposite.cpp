#ifndef TEXT_PARAM_COMPOSITE_CPP
#define TEXT_PARAM_COMPOSITE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "TextComposite.h"   // header
#include <util/util/ioUtil.h>    

#include <vector>

namespace Util
{
  
   /** 
   * Constructor
   */
   TextComposite::TextComposite()
   {}
  
   /** 
   * Virtual destructor.
   */
   TextComposite::~TextComposite()
   {}

   /**
   * Resets TextComposite to its empty state.
   */
   void TextComposite::resetParam()
   {  lines_.clear(); }

   /** 
   * Read all parameters from an input stream.
   */
   void TextComposite::readParam(std::istream &in)
   {
      std::string line;
      int len;
      int depth = 0;

      getline(in, line);
      len = rStrip(line);
      if (line[len-1] == '{') {
         ++depth;
      } else {
         UTIL_THROW("Malformed beginning line");
      }
      lines_.push_back(line);

      while (depth > 0) {
         getline(in, line);
         len = rStrip(line);
         lines_.push_back(line);
         if (len > 0) {
            if (line[len-1] == '{') {
               ++depth;
            }
            if (line[len-1] == '}') {
               --depth;
            }
         }
      }
    
      /// Find the number of leading spaces 
      size_t found;
      size_t min = 1028;
      char space = ' ';
      for (size_t i = 0; i < lines_.size(); ++i) {
         found = lines_[i].find_first_not_of(space); 
         if (found < min) {
            min = found;
         }
      }

      // Strip min leading spaces from each line
      if (min > 0) {
         for (size_t i = 0; i < lines_.size(); ++i) {
            lines_[i].erase(0, min);
         }
      }

   }

   /*
   * Write all parameters to an output stream.
   */
   void TextComposite::writeParam(std::ostream &out)
   {
      int size = lines_.size();
      for (int i = 0; i < size; ++i) {
         out << indent() << lines_[i] << std::endl;
      }
   }
   
}
#endif
