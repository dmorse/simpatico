#ifndef MCMD_COMMAND_H
#define MCMD_COMMAND_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>
#include <util/random/Random.h>
#include <util/global.h>

namespace McMd
{

   using namespace Util;

   class Simulation;

   /**
   * Command is an object that can be invoked from the command script.
   *
   * Usage: 
   *
   * Usage: In this snippet, "in" is a std::istream with a cursor set at the beginning
   * of a line in command script, and "command" is an instance of a subclass of the 
   * Command base classe, representing a specific command.
   * \code
   *  
   *   std::string name;
   *   in >> name;
   *   if (command.match(name)) {
   *      command.execute(in);
   *   }
   *
   * \endcode
   * The readCommand() function of a CommandManager instead applies the match command
   * to a list of commands and executes the first one with a name that matches.
   *
   * \ingroup McMd_Command_Module
   */
   class Command : public ParamComposite
   {

   public:

      /**
      * Constructor.
      *
      * \param name identifier for this command, as the first word in a command line.
      */
      Command(std::string name);

      /**
      * Destructor.
      *
      * Empty default implementation.
      */
      virtual ~Command();

      /**
      * Compare a string to the name of this command.
      *
      * \return true if name matches the name for this command, false otherwise.
      */
      virtual bool match(std::string name);

      /**
      * Execute this command.
      *
      * This function should be called after the command name is read from the
      * istream in, and the match() function has been called and returned true.
      * The std::itream in should be the input stream from which the name was 
      * read, with the cursor initially set after the end of the name string.
      *
      * \param in input stream from which the command name was read.
      */
      virtual void execute(std::istream& in) = 0;

      /**
      * Output statistics accumulated by this command (called at the end of the simulation).
      *
      * Empty default implementation.
      */
      virtual void output();

   private:

      /// Identifier for this command in a command script
      std::string name_;

   };

}
#endif
