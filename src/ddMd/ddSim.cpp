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

   bool  eflag  = false;
   char* rarg   = 0;

   // Read command-line arguments
   int c;
   opterr = 0;
   while ((c = getopt(argc, argv, "epr:")) != -1) {
      switch (c) {
      case 'e':
        eflag = true;
        break;
      case '?':
        std::cout << "Unknown option -" << optopt << std::endl;
        return 1;
      }
   }

   #ifdef UTIL_MPI
   MPI::Init();
   DdMd::Simulation simulation(MPI::COMM_WORLD);
   #else
   DdMd::Simulation simulation();
   #endif

   if (eflag) {
      // Enable echoing of parameters to log file as they are read.
      Util::ParamComponent::setEcho(true);
   }

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
