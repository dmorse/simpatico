/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/simulation/Simulation.h>
#include <util/param/ParamComponent.h>
#include <util/global.h>

#include <iostream>

using namespace Util;

/**
* Program for parallel domain-decomposition molecular dynamics simulation.
*
* Command line options:
*
*   -e  
*    Enable echoing of parameter file to log file as it is read. 
*    This is sometimes useful for debugging the parameter file.
*
*  -s nSystem 
*   Split communicator nSystem processor groups, each for a different system.
*/
int main(int argc, char **argv)
{

   #ifdef UTIL_MPI
   MPI::Init();
   #endif

   #ifdef UTIL_MPI
   DdMd::Simulation simulation(MPI::COMM_WORLD);
   #else
   DdMd::Simulation simulation();
   #endif

   // Read and apply command line options.
   simulation.setOptions(argc, argv); 

   // Read parameter file from default param stream (read on master).
   simulation.readParam(); 

   // Read command file (path is specified in parameter file).
   simulation.readCommands();

   #ifdef UTIL_MPI
   MPI::Finalize();
   #endif

   return 0;
}
