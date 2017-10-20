#ifndef MCMD_MD_DEFORM_COMMAND_H
#define MCMD_MD_DEFORM_COMMAND_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/commands/DeformCommand.h>
#include <mcMd/mdSimulation/MdSystem.h>

namespace McMd
{

   class MdSystem;
   class MdPairPotential;
   using namespace Util;

   /**
   * Command to deform the unit cell.
   *
   * \ingroup McMd_Command_Module
   */
   class MdDeformCommand : public DeformCommand
   {

   public:

      /**
      * Constructor.
      */
      MdDeformCommand(MdSystem& system);

      /**
      * Destructor.
      */
      virtual ~MdDeformCommand();

      /**
      * Rebuild cell list.
      */
      virtual void reneighbor();

   private:

      MdPairPotential* pairPtr_;

   };

}
#endif
