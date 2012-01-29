#ifndef OUT_FILE_H
#define OUT_FILE_H

#include <fstream>

namespace Util
{

   /**
   * A std::ofstream that stores a filename.
   *
   */
   class OutFile : public std::ofstream
   {
   
   public:

      /**
      * Default constructor.
      */   
      OutFile();
   
      /**
      * Construct and open file, and store filename.
      */   
      explicit OutFile(const char* filename, 
                       ios_base::openmode mode = ios_base::out);
   
      /**
      * Open the file, and store the filename.
      */   
      void open(const char* filename);
   
      /**
      * Close the file.
      */
      void close();
   
      /**
      * Access the file name by reference.
      *
      * Because this function returns by non-const reference, its value may
      * be used at a left-value to modify the filename.
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
