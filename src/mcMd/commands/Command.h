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
   * \ingroup McMd_Command_Module
   */
   class Command : public ParamComposite
   {

   public:

      /**
      * Constructor.
      */
      Command();

      /**
      * Destructor.
      *
      * Empty default implementation.
      */
      virtual ~Command();

      /**
      * Read and execute command.
      * 
      * The string name is the capitalized command name that begins a command.
      * If the name string matches the command name for this command, the rest 
      * of the line is read, the command is executed, and the function returns
      * true. If the name does not match, the function immediately returns 
      * false without reading from istream in.
      *
      * \return true if name is recognized, false if not
      */
      virtual bool readCommand(std::string const & name, std::istream& in) = 0;

      /**
      * Output statistics for this command (called at the end of the simulation)
      */
      virtual void output();

   };

}
#endif
