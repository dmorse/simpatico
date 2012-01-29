#ifndef OUT_FILE_CPP
#define OUT_FILE_CPP

#include "Outfile.h"

namespace Util
{

   /// Default constructor.
   OutFile::OutFile()
    : std::ofstream(),
      filename_()
   {}

   /// Constructor that opens file.
   OutFile::OutFile(const char* filename, ios_base::openmode mode)
    : std::ofstream(filename, mode),
      filename_(filename)
   {}

   /// Open file.
   void OutFile::open(const char* filename)
   {
      std::ofstream::open(filename);
      filename_ = filename;
   }

   /// Close file.
   void OutFile::close()
   {
      std::ofstream::close();
      filename_ = "";
   }

   /// Get filename.
   const std::string& OutFile::filename() const
   {  return filename_; }
   
}
#endif
