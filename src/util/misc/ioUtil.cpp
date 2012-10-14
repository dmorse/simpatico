#ifndef UTIL_IO_UTIL_CPP
#define UTIL_IO_UTIL_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "ioUtil.h"
#include "Log.h"
#include <util/global.h>

namespace Util
{

   /*
   * Strip trailing whitespace from a string.
   */
   int rStrip(std::string& str)
   {
      size_t found;
      std::string whitespaces(" \t\n\r");
      found = str.find_last_not_of(whitespaces);
      if (found != std::string::npos) {
        str.erase(found + 1);
        return int(found + 1);
      } else {
        str.clear();
        return 0;
      }
   }

   /*
   * Read string, and compare to expected value. 
   *
   * Throw Exception if input value differs from expected value.
   */
   void checkString(std::istream& in, std::string& expected)
   {
      std::string actual;
      in >> actual;
      if ( actual != expected ) {
         Log::file() << "Error in checkString"     << std::endl;
         Log::file() << "Expected: " <<  expected  << std::endl;
	 Log::file() << "Scanned:  " <<  actual    << std::endl;
         UTIL_THROW("Incorrect string");
      };
   }

   /*
   * Return std::string representation of an integer.
   */
   std::string toString(int n)
   {
      std::stringstream ss;
      ss << n;
      return ss.str();
   }


   /*
   * Get the next line of input skipping an blank lines.
   */
   void getNextLine(std::istream& in, std::string& line)
   {
      bool hasLine = false;
      while (!hasLine) {
         std::getline(in, line);
         if (!line.empty()) {
            hasLine = true;
         }
      }
   }

   /*
   * Transfer the next line of input to a stringstream
   */
   void getLine(std::istream& in, std::stringstream& line)
   {
      std::string string;
      line.str(string);
      getline(in, string);
      line.str(string);
   }

   /*
   * Transfer the next line of input to an sstream
   */
   void getNextLine(std::istream& in, std::stringstream& line)
   {
      std::string string;
      line.str(string);
      bool hasLine = false;
      while (!hasLine) {
         std::getline(in, string);
         if (!string.empty()) {
            hasLine = true;
         }
      }
      line.str(string);
   }

}
#endif
