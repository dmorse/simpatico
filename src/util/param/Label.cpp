#ifndef UTIL_LABEL_CPP
#define UTIL_LABEL_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Label.h"
#include <util/global.h>

#include <iomanip>

//using namespace std;

namespace Util
{

   /*
   * Constructor.
   */
   Label::Label(const char* label)
    : label_(label)
   {}

   /*
   * Copy constructor.
   */
   Label::Label(const Label& other)
    : label_(other.label_)
   {}

   /*
   * Destructor.
   */
   Label::~Label()
   {}

   /*
   * Return label string.
   */
   std::string Label::string() const
   {  return label_; }

   /*
   * Extract a label from an input stream.
   */
   std::istream& operator>>(std::istream& in, Label label)
   {
      std::string actual;
      in >> actual;
      if ( actual != label.label_ ) {
         Log::file() << "Error reading label"          << std::endl;
         Log::file() << "Expected: " <<  label.label_  << std::endl;
	 Log::file() << "Scanned:  " <<  actual        << std::endl;
         UTIL_THROW("Incorrect label");
      };
      return in;
   }

   /*
   * Insert a Label into an output stream.
   */
   std::ostream& operator<<(std::ostream& out, Label label)
   {
      out << std::left << std::setw(Label::LabelWidth) 
          << label.label_; 
      out << std::right;
      return out;
   }

} 
#endif
