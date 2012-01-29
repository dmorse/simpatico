#include <mcMd/mdSimulation/MdSimulation.h>
#include <util/param/ParamComponent.h>

#include <unistd.h>
#include <memory>

using namespace McMd;
using namespace Util;

/**
* Main program for Monte Carlo Simulation.
*
* Options:
*
*   -e  Enable echoing of the parameter file to the log file as it
*       is read.
*
*   -p  Enable use of a free energy perturbation. 
*
*   -r filename
*      Restart a simulation. The required parameter "filename" is 
*      the base name for 3 input files: filename.prm, filename.rst, 
*      and filename.cmd.
*
* Files:
*
* If compiled in serial mode (ifndef UTIL_MPI), the default parameter 
* stream read by readParam() is std::cin and the default log file 
* Log::file() is std::cout.
*
* If compiled in parallel mode (ifdef UTIL_MPI) and invoked with the 
* "-p" option, a single parameter file is read from std::cin, as in a 
* serial job.  If compiled in parallel mode and invoked without the -p 
* option, each processor reads from a different parameter file, which 
* is named "n/param" for processor n.
*
* If compiled in parallel mode (ifdef UTIL_MPI), paths for all output
* files produced by the processor with rank n are prefixed by "n/". All
* output files produced by processor 6, for example, are thus put in a 
* directory named "6/".  The default log file for processor n is "n/log". 
*/
int main(int argc, char **argv)
{

   bool  eflag  = false;
   bool  rflag  = false;
   char* rarg   = 0;
   #ifdef MCMD_PERTURB
   bool pflag = false;
   #endif

   // Read program arguments
   int c;
   opterr = 0;
   while ((c = getopt(argc, argv, "epr:")) != -1) {
      switch (c) {
      case 'e':
        eflag = true;
        break;
      case 'r':
        rflag = true;
        rarg  = optarg;
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

   MdSimulation simulation;

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
      simulation.setParamCommunicator();
      #endif

   }
   #endif

   if (rflag) {

      std::cout << "Reading restart" << std::endl;
      std::cout << "Base file name " << std::string(rarg) << std::endl;

      simulation.readRestart(std::string(rarg));

   } else {

      // Read parameters from default parameter file
      simulation.readParam();
   
      // Read command script to run simulation
      simulation.readCommands();

   }

   // Normal completion
   return 0;

}
