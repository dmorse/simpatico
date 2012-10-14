#ifndef UTIL_BLANK_CPP
#define UTIL_BLANK_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Blank.h"

namespace Util
{

   /// Constructor.
   Blank::Blank()
    : ParamComponent()
   {}

   /// Virtual Destructor
   Blank::~Blank()
   {}

   /**
   * Read a blank line
   *
   * \param in input stream
   */
   void Blank::readParam(std::istream &in)
   {
      if (isParamIoProcessor()) {
         char buf[255];
         in.getline(buf,255);
         if (ParamComponent::echo()) {
            writeParam(Log::file());
         }
      }
   }

   /**
   * Write a blank line
   *
   * \param out output stream
   */
   void Blank::writeParam(std::ostream &out)
   {  out << std::endl; }

} 
#endif
