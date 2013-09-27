#ifndef UTIL_BEGIN_CPP
#define UTIL_BEGIN_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Begin.h"
#include "Label.h"
#include <util/global.h>

namespace Util
{

   /* 
   * Begin constructor
   */
   Begin::Begin(const char *label)
    : label_(label)
   {} 

   /* 
   * Read label, check against expected value
   */
   void Begin::readParam(std::istream &in)
   {
      if (isIoProcessor()) {
         std::string expected;
         expected = label_ + "{";
         in >> Label(expected.c_str());
         if (ParamComponent::echo()) {
            writeParam(Log::file());
         }
      }
   }

   /* 
   * Begin::writeParam() template
   */
   void Begin::writeParam(std::ostream &out)
   {  out << indent() << label_ << "{" << std::endl; }

   /*
   * Do-nothing implementation of virtual resetIo function.
   */
   void Begin::resetParam()
   {}

} 
#endif
