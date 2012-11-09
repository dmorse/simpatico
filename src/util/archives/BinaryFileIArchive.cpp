#ifndef UTIL_BINARY_FILE_I_ARCHIVE_CPP
#define UTIL_BINARY_FILE_I_ARCHIVE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "BinaryFileIArchive.h"

#include <vector>

namespace Util
{

   /*
   * Constructor.
   */
   BinaryFileIArchive::BinaryFileIArchive()
    : istreamPtr_(0),
      version_(0)
   {}

   /*
   * Destructor.
   */
   BinaryFileIArchive::~BinaryFileIArchive()
   {}  

   /*
   * Allocate a block of memory.
   */
   void BinaryFileIArchive::setStream(std::istream& out)
   { 
      istreamPtr_ = &out;
   }

   /*
   * Load a std::string from BinaryFileIArchive.
   */
   template <>
   void serialize(BinaryFileIArchive& ar, std::string& data, 
                  const unsigned int version)
   {
      static std::vector<char> charvec;
      size_t size;
      ar.unpack(size);
      if (size > charvec.capacity()) {
         charvec.reserve(size + 8);
      }
      ar.unpack(&charvec[0], size);
      data = &charvec[0];
   }


}
#endif
