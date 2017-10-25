/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "MdDeformCommand.h"
#include <mcMd/mdSimulation/MdSystem.h>
#include <mcMd/potentials/pair/MdPairPotential.h>

namespace McMd
{

   using namespace Util;

   /*
   * Default constructor.
   */
   MdDeformCommand::MdDeformCommand(MdSystem& system) 
    : DeformCommand(system),
      pairPtr_(&system.pairPotential())
   {  setClassName("MdDeformCommand"); }

   /*
   * Default destructor.
   */
   MdDeformCommand::~MdDeformCommand()
   {}

   /*
   * Rebuild pair list after deformation.
   */
   void MdDeformCommand::reneighbor()
   {
      #ifndef SIMP_NOPAIR 
      pairPtr_->buildPairList();
      #endif
   }

}
