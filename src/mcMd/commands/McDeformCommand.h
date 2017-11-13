#ifndef MCMD_MC_DEFORM_COMMAND_H
#define MCMD_MC_DEFORM_COMMAND_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/commands/DeformCommand.h>
#include <mcMd/mcSimulation/McSystem.h>

namespace McMd
{

   class McSystem;
   class McPairPotential;
   using namespace Util;

   /**
   * Command to deform the unit cell.
   *
   * \ingroup McMd_Command_Module
   */
   class McDeformCommand : public DeformCommand
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
      * Rebuild cell list.
      */
      virtual void reneighbor();

   private:

      #ifndef SIMP_NOPAIR
      McPairPotential* pairPtr_;
      #endif

   };

}
#endif
