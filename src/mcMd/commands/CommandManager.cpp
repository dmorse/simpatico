/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "CommandManager.h"                

namespace McMd
{

   using namespace Util;

   // Constructor.
   CommandManager::CommandManager()
    : Manager<Command>(true)
   {  setClassName("CommandManager"); }


   // Destructor.
   CommandManager::~CommandManager()
   {}

   bool 
   CommandManager::readCommand(std::string name, std::istream& in) 
   {
      // Attempt to match one of the standard commands
      bool success = readStandardCommand(name, in);
      if (success) {
         return true;
      }

      // Attempt to match a customized command
      for (int i = 0; i < size(); ++i) {
         Command& command = (*this)[i];
         if (command.match(name)) {
            command.execute(in);
            // Success: Match and execution
            return true;
         }
      }

      // Failure: If this point was reach, no match was found
      return false;
   }

   bool 
   CommandManager::readStandardCommand(std::string name, std::istream& in) 
   { return false; }

}
