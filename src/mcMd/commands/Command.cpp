/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Command.h"
#include <util/archives/Serializable_includes.h>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   Command::Command(std::string name) 
    : name_(name)
   {}

   /*
   * Default destructor.
   */
   Command::~Command()
   {}

   /*
   * Attempt to match a string to the name identifier for this command.
   */
   bool Command::match(std::string name)
   {  return (name == name_); }

}
