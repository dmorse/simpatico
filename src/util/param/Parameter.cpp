#ifndef UTIL_PARAMETER_CPP
#define UTIL_PARAMETER_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Parameter.h"

namespace Util
{

   /*
   * Constructor. 
   */
   Parameter::Parameter(const char *label)
    : label_(label)
   {}

   /*
   * Destructor.
   */
   Parameter::~Parameter()
   {}

   /*
   * Return label string.
   */
   std::string Parameter::label() 
   {  return label_.string(); }


} 
#endif
