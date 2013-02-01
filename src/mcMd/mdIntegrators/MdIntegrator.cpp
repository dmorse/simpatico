#ifndef MCMD_MD_INTEGRATOR_CPP
#define MCMD_MD_INTEGRATOR_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "MdIntegrator.h"
#include <mcMd/mdSimulation/MdSystem.h>
#include <mcMd/simulation/Simulation.h>
#include <util/archives/Serializable_includes.h>


namespace McMd
{

   using namespace Util;

   /* 
   * Constructor.   
   */
   MdIntegrator::MdIntegrator(MdSystem& system)
   : dt_(0.0),
     boundaryPtr_(&system.boundary()),
     systemPtr_(&system),
     simulationPtr_(&system.simulation())
   {}

   /* 
   * Destructor.   
   */
   MdIntegrator::~MdIntegrator() 
   {}

}
#endif
