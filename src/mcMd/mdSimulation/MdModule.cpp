#ifndef MCMD_MD_MODULE_CPP
#define MCMD_MD_MODULE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "MdModule.h"
#include <mcMd/mdSimulation/MdSimulation.h>

#include <util/global.h>

namespace McMd
{

   /*
   * Constructor.
   */
   MdModule::MdModule(MdSimulation& simulation) 
    : simulationPtr_(&simulation),
      systemPtr_(&simulation.system())
   {}

   /*
   * Destructor.
   */
   MdModule::~MdModule()
   {}

   /*
   * Return MdSimulation by reference.
   */
   MdSimulation& MdModule::simulation() const
   {   return *simulationPtr_; }

   /*
   * Return MdSystem by reference.
   */
   MdSystem& MdModule::system() const
   {   return *systemPtr_; }

}
#endif
