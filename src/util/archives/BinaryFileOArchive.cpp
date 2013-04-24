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
    : filePtr_(0),
      version_(0),
      createdFile_(true)
   {  filePtr_ = new std::ofstream(); }

   /*
   * Constructor.
   */
   BinaryFileOArchive::BinaryFileOArchive(std::string filename)
    : filePtr_(0),
      version_(0),
      createdFile_(true)
   {  filePtr_ = new std::ofstream(filename.c_str()); }

   /*
   * Constructor.
   */
   BinaryFileOArchive::BinaryFileOArchive(std::ofstream& file)
    : filePtr_(&file),
      version_(0),
      createdFile_(false)
   {}

   /*
   * Destructor.
   */
   BinaryFileOArchive::~BinaryFileOArchive()
   {
      if (filePtr_ && createdFile_) {  
         delete filePtr_; 
      }
   }  

   /*
   * Return underlying file by reference.
   */
   std::ofstream& BinaryFileOArchive::file()
   {  return *filePtr_; }

}
#endif
