#ifndef MCMD_DEFORM_COMMAND_H
#define MCMD_DEFORM_COMMAND_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/commands/Command.h>
#include <mcMd/mcSimulation/McSystemInterface.h>
#include <mcMd/mcSimulation/McSystem.h>

namespace McMd
{

   using namespace Util;

   /**
   * Command to deform the unit cell.
   *
   * \ingroup McMd_Command_Module
   */
   class McDeformCommand : public Command, private McSystemInterface
   {

   public:

      /**
      * Constructor.
      */
      McDeformCommand(McSystem& system);

      /**
      * Destructor.
      */
      virtual ~McDeformCommand();

      /**
      * Read and execute command iff name == "DEFORM_CELL".
      *
      * \return true if name is recognized, false if not
      */
      virtual 
      bool readCommand(std::string const & name, std::istream& in);

   };

}
#endif
