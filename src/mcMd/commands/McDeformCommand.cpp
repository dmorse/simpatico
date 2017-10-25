/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McDeformCommand.h"
#include <mcMd/mcSimulation/McSystem.h>
#include <mcMd/potentials/pair/McPairPotential.h>

namespace McMd
{

   using namespace Util;

   /*
   * Default constructor.
   */
   McDeformCommand::McDeformCommand(McSystem& system) 
    : DeformCommand(system),
      pairPtr_(&system.pairPotential())
   {  setClassName("McDeformCommand"); }

   /*
   * Default destructor.
   */
   McDeformCommand::~McDeformCommand()
   {}

   /*
   * Rebuild cell list after deformation.
   */
   void McDeformCommand::reneighbor()
   {
      #ifndef SIMP_NOPAIR 
      // Generate cell list
      pairPtr_->buildCellList();
      #endif
   }

}
