#ifndef MCMD_DEFORM_COMMAND_H
#define MCMD_DEFORM_COMMAND_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/commands/Command.h>
#include <mcMd/simulation/System.h>
#include <mcMd/simulation/SystemInterface.h>

namespace McMd
{

   using namespace Util;

   /**
   * Command to deform the unit cell.
   *
   * \ingroup McMd_Command_Module
   */
   class DeformCommand : public Command, private SystemInterface
   {

   public:

      /**
      * Constructor.
      */
      DeformCommand(System& system);

      /**
      * Destructor.
      */
      virtual ~DeformCommand();

      /**
      * Execute deformation command.
      */
      virtual 
      void execute(std::istream& in);

   protected:

      /** 
      * Rebuild cell and/or pair list after deforming positions.
      */
      virtual void reneighbor() = 0;

   };

}
#endif
