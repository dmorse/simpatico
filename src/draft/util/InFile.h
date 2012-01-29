#ifndef IN_FILE_H
#define IN_FILE_H

#include <fstream>

namespace Util
{

   /**
   * An improved std::ifstream.
   *
   * An Infile is an ifstream that stores its own filename,
   * and that automatically handles errors encountered when
   * a file is opened. 
   */
   class InFile : public std::ifstream
   {
   
   public:

      /**
      * Default constructor.
      */   
      InFile();
   
      /**
      * Construct and open file. 
      */   
      explicit InFile(const char* filename, 
                       ios_base::openmode mode = ios_base::in);
   
      /**
      * Open the file.
      */   
      void open(const char* filename);
   
      /**
      * Close the file.
      */
      void close();
   
      /**
      * Get the file name.
      *
      * Return the filename used to open this file if it is open, or an 
      * empty string if it is not open.
      */   
      const std::string& filename() const;

   private:
  
      /**
      * Name of file if the file is open, or an empty string otherwise.
      */ 
      std::string filename_;
   
   };
   
}
#endif
