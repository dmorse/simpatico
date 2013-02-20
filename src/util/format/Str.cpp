#ifndef UTIL_STR_CPP
#define UTIL_STR_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Str.h"
#include "Format.h"

#include <iostream>

namespace Util
{

   /// Default constructor.
   Str::Str()
    : value_(),
      width_(Format::defaultWidth())
   {}

   /// Constructor, value only (explicit).
   Str::Str(std::string value)
    : value_(value),
      width_(Format::defaultWidth())
   {}

   /// Constructor, value and width.
   Str::Str(std::string value, int width)
    : value_(value),
      width_(width)
   {}

   void Str::setValue(std::string value)
   {  value_ = value; }

   void Str::setWidth(int width)
   {  width_ = width; }

   std::string Str::value() const
   {  return value_; }

   int Str::width() const
   {  return width_; }

   /*
   * Input stream extractor for an Str object.
   */
   std::istream& operator>>(std::istream& in, Str &object)
   {
      in >> object.value_;
      return in;
   }

   /*
   * Output stream inserter for an Str object.
   */
   std::ostream& operator<<(std::ostream& out, const Str &object)
   {
      out.width(object.width_);
      out << object.value_;
      return out;
   }

} 
#endif
