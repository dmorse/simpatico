#ifndef MCMD_COMMAND_MANAGER_H
#define MCMD_COMMAND_MANAGER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
w Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Manager.h>      // base class template
#include <mcMd/commands/Command.h>   // class template argument

namespace McMd
{

   using namespace Util;

   /**
   * Manager for Command objects in an MdSimulation.
   *
   * \ingroup McMd_Manager_Module
   * \ingroup McMd_Command_Module
   */
   class CommandManager : public Manager<Command>
   {

   public:

      /**
      * Constructor.
      */
      CommandManager();

      /**
      * Destructor.
      */
      virtual ~CommandManager();

      /**
      * Attempt to read a command line, execute if recognized.
      *
      * \return true iff succesful match to name
      */
      bool readCommand(std::string const & name, std::istream& in);

   };

}
#endif
