#include "../../src/mcSimulation/McSimulation.h"
#include "../../src/mcSimulation/McSystem.h"
#include "../../src/simulation/ConfigIo.h"

#include "../../src/configIos/PmcConfigIo.h"
#include "../../src/configIos/LammpsConfigIo.h"

using namespace McMd;
using namespace Util;

#include <cstdio>
#include <iostream>
#include <ctime>

/**
* Example of a program for converting a config file format.
*
* This program instantiates and uses an McSimulation object.
*
* \ingroup Example_Module
*/
int main()
{

   McSimulation    simulation;

   // Create and register a PmcConfigIo
   // PmcConfigIo pmcConfigIo(simulation.system());
   // simulation.system().setConfigIo(pmcConfigIo);

   // Read parameter and config files 
   simulation.readParam(std::cin);
   simulation.writeParam(std::cout);

   // Create and register a default ConfigIo
   ConfigIo configIo(simulation.system());
   simulation.system().setConfigIo(configIo);

   // Create and register a LammpsConfigIo
   LammpsConfigIo lammpsConfigIo(simulation.system());
   simulation.system().setConfigIo(lammpsConfigIo);

   // Write config file
   simulation.system().writeConfig();

}
