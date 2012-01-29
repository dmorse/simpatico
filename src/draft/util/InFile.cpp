#ifndef IN_FILE_CPP
#define IN_FILE_CPP

#include "Infile.h"

namespace Util
{

   /// Default constructor.
   InFile::InFile()
    : std::ifstream(),
      filename_()
   {}

   /// Constructor that opens file.
   InFile::InFile(const char* filename, ios_base::openmode mode)
    : std::ifstream(filename, mode),
      filename_(filename)
   {}

   /// Open file.
   void InFile::open(const char* filename)
   {
      std::ifstream::open(filename);
      filename_ = filename;
   }

   /// Close file.
   void InFile::close()
   {
      std::ifstream::close();
      filename_ = "";
   }

   /// Get filename.
   const std::string& InFile::filename() const
   {  return filename_; }
   
}
#endif
