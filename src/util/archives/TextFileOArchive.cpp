#ifndef UTIL_TEXT_FILE_O_ARCHIVE_CPP
#define UTIL_TEXT_FILE_O_ARCHIVE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "TextFileOArchive.h"

namespace Util
{

   /*
   * Constructor.
   */
   TextFileOArchive::TextFileOArchive()
    : ostreamPtr_(0),
      version_(0)
   {}

   /*
   * Destructor.
   */
   TextFileOArchive::~TextFileOArchive()
   {}  

   /*
   * Allocate a block of memory.
   */
   void TextFileOArchive::setStream(std::ostream& out)
   { 
      ostreamPtr_ = &out;
   }

}
#endif
