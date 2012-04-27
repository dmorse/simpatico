#ifndef DDMD_DDSIM_CPP
#define DDMD_DDSIM_CPP

#include <ddMd/simulation/Simulation.h>
#include <util/param/ParamComponent.h>

/**
* Main program for paralel domain-decomposition MD simulations.
*
* Options:
*
*   -e  Enable echoing of the parameter file to the log file as it
*       is read.
*/
int main(int argc, char **argv)
{

   #ifdef UTIL_MPI
   MPI::Init();
   DdMd::Simulation simulation(MPI::COMM_WORLD);
   #else
   DdMd::Simulation simulation();
   #endif

   // Read parameter file from standard input (read on master).
   simulation.setOptions(argc, argv); 

   // Read parameter file from standard input (read on master).
   simulation.readParam(std::cin); 

   // Read command file (path is specified in parameter file).
   simulation.readCommands();

   #ifdef UTIL_MPI
   MPI::Finalize();
   #endif

   return 0; // normal completion
}
#endif
