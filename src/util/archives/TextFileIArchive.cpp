#ifndef UTIL_TEXT_FILE_I_ARCHIVE_CPP
#define UTIL_TEXT_FILE_I_ARCHIVE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "TextFileIArchive.h"

#include <vector>

namespace Util
{

   /*
   * Constructor.
   */
   TextFileIArchive::TextFileIArchive()
    : istreamPtr_(0),
      version_(0)
   {}

   /*
   * Destructor.
   */
   TextFileIArchive::~TextFileIArchive()
   {}  

   /*
   * Allocate a block of memory.
   */
   void TextFileIArchive::setStream(std::istream& out)
   { 
      istreamPtr_ = &out;
   }

   /*
   * Load a std::string from TextFileIArchive.
   */
   template <>
   void serialize(TextFileIArchive& ar, std::string& data, 
                  const unsigned int version)
   {
      static std::vector<char> charvec;
      size_t size;
      ar.unpack(size);
      if (size > charvec.capacity()) {
         charvec.reserve(size + 8);
      }
      // Read endline character after size
      char   endline;
      ar.unpack(endline);
      // Read actual string.
      ar.unpack(&charvec[0], size);
      data = &charvec[0];
   }


}
#endif
