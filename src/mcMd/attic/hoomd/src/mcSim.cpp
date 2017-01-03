#include <mcMd/mcSimulation/McSimulation.h>
#include <util/param/ParamComponent.h>

// module definition
#include "HoomdMcModule.h"

// Custom factory classes
#include "mcMoves/HoomdMcMoveFactory.h"
#include "potentials/pair/HoomdPairFactory.h"
#include "potentials/external/HoomdExternalFactory.h"

#include <unistd.h>
#include <memory>

#ifdef UTIL_MPI
#include <mpi.h>
#endif

using namespace McMd;
using namespace Util;

/**
* Main program for Monte Carlo Simulation.
*
* Options:
*
*   -e  Enable echoing of the parameter file to the log file while
*       the parameter file is read. This can be helpful if an error
*       is encountered while the parameter file is being read, to
*       help locate the syntax error in the parameter file.
*
*   -p  Enable use of a thermodynamic perturbation. A perturbation
*       defines how a Hamiltonian or (more generally) a Boltzmann
*       weight depends on a parameter. If used with a program that
*       is compiled for parallel operation (ifdef UTIL_MPI), this
*       option causes a master processor to read a single parameter
*       file, broadcast a set of default parameters to all of the 
*       processors, and then assign each processor a different value 
*       of a perturbation parameter that is used to modify the state 
*       on each processor. The simplest case is when the perturbation 
*       parameter may be the inverse temperature (beta = 1/kT), in
*       which case different processors simulate similar system at 
*       different temperatures.
*
* Files:
*
* If compiled in serial mode (ifndef UTIL_MPI), the default parameter 
* stream read by readParam() is std::cin and the default log file 
* Log::file() is std::cout.
*
* If compiled in parallel mode (ifdef UTIL_MPI), the paths for all output
* files produced by the processor with rank n are prefixed by "n/". All
* output files produced by processor 6, for example, are thus put in a 
* directory named "6/".  The default log file for processor n is "n/log". 
* If compiled in parallel mode and invoked with the "-p" option, to enable 
* a perturbation, the default parameter stream is std::cin. In this case, 
* this file is read by the master processor and broadcast to other 
* processors. If compiled in parallel mode and invoked without the -p 
* option, each processor reads from a different parameter file, which is 
* named "n/param" for processor n.
*/
int main(int argc, char **argv)
{

   bool eflag = false;
   #ifdef MCMD_PERTURB
   bool pflag = false;
   #endif

   // Read program arguments
   int c;
   opterr = 0;
   while ((c = getopt(argc, argv, "pe")) != -1) {
      switch (c) {
      case 'e':
        eflag = true;
        break;
      #ifdef MCMD_PERTURB
      case 'p':
        pflag = true;
        break;
      #endif
      case '?':
        std::cout << "Unknown option -" << optopt << std::endl;
        return 1;
      }
   }

   #ifdef UTIL_MPI
   MPI::Init();
   McMd::McSimulation simulation(MPI::COMM_WORLD);
   #else
   McMd::McSimulation simulation;
   #endif

   // Set flag to echo parameters as they are read.
   if (eflag) {
      ParamComponent::setEcho(true);
   }

   #ifdef MCMD_PERTURB

   // Set to use a perturbation.
   if (pflag) {

      // Set to expect perturbation in the param file.
      simulation.system().setExpectPerturbation();

      #ifdef UTIL_MPI
      Log::file() << "Set to read parameters from a single file" << std::endl;
      simulation.setIoCommunicator();
      #endif

   }
   #endif

   // set up module
   HoomdMcModule module(simulation);
   module.addFactories();
 
   // Read parameters from default parameter file
   simulation.readParam();
 
   // Read command script to run simulation
   simulation.readCommands();

   #ifdef UTIL_MPI
   MPI::Finalize();
   #endif
   // Normal completion
   return 0;

}
