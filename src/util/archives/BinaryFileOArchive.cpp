#ifndef UTIL_BINARY_FILE_O_ARCHIVE_CPP
#define UTIL_BINARY_FILE_O_ARCHIVE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "BinaryFileOArchive.h"

namespace Util
{

   /*
   * Constructor.
   */
   BinaryFileOArchive::BinaryFileOArchive()
    : ostreamPtr_(0),
      version_(0)
   {}

   /*
   * Destructor.
   */
   BinaryFileOArchive::~BinaryFileOArchive()
   {}  

   /*
   * Allocate a block of memory.
   */
   void BinaryFileOArchive::setStream(std::ostream& out)
   { 
      ostreamPtr_ = &out;
   }

}
#endif
