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

namespace Util
{

   bool  Label::isMatch_ = true;
   std::string  Label::input_;

   /*
   * Clear read buffer (static variable).
   */
   void Label::clear()
   {
      isMatch_ = true;
      input_.clear(); 
   }

   /*
   * Constructor.
   */
   Label::Label(const char* label, bool isRequired)
    : isRequired_(isRequired),
      label_(label)
   {}

   /*
   * Copy constructor.
   */
   Label::Label(const Label& other)
    : isRequired_(other.isRequired_),
      label_(other.label_)
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
      // If previous input value matched, read a new one.
      if (label.isMatch_) {
         in >> label.input_;
      }
      if (label.input_ == label.label_) {
         label.isMatch_ = true;
         label.input_.clear();
      } else {
         label.isMatch_ = false;
         if (label.isRequired_) {
            Log::file() << "Error reading label"        << std::endl;
            Log::file() << "Expected: " << label.label_ << std::endl;
            Log::file() << "Scanned:  " << label.input_ << std::endl;
            UTIL_THROW("Incorrect label");
         }
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
