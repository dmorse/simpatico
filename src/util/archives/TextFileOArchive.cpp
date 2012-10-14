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
    : filePtr_(0),
      version_(0)
   {  filePtr_ = new std::ofstream(); }

   /*
   * Constructor.
   */
   TextFileOArchive::TextFileOArchive(std::string filename)
    : filePtr_(0),
      version_(0)
   {  filePtr_ = new std::ofstream(filename.c_str()); }

   /*
   * Destructor.
   */
   TextFileOArchive::~TextFileOArchive()
   {  delete filePtr_; }  

   /*
   * Return underlying file by reference.
   */
   std::ofstream& TextFileOArchive::file()
   {  return *filePtr_; }

}
#endif
