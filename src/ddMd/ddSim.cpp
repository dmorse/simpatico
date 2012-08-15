/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/simulation/Simulation.h>
#include <util/param/ParamComponent.h>

/**
* Program for parallel domain-decomposition molecular dynamics simulation.
*
* Command line options:
*
*   -e  
*    Enable echoing of parameter file to log file as it is read.
*/
int main(int argc, char **argv)
{

   #ifdef UTIL_MPI
   MPI::Init();
   DdMd::Simulation simulation(MPI::COMM_WORLD);
   #else
   DdMd::Simulation simulation();
   #endif

   // Read command line options (read on master).
   simulation.setOptions(argc, argv); 

   // Read parameter file from standard input (read on master).
   simulation.readParam(std::cin); 

   // Read command file (path is specified in parameter file).
   simulation.readCommands();

   #ifdef UTIL_MPI
   MPI::Finalize();
   #endif

   return 0;
}
