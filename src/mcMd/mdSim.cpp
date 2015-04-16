/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/mdSimulation/MdSimulation.h>

/**
* \page mdSim_page mdSim - serial molecular dynamics program
*
* Single-processor molecular dynamics simulation program.
*
* Usage:
*
*    mdSim [-e] [-r file] [-p file] [-c file] [-i prefix] [-o prefix] [-f]
*
* Options:
*
*   -e  
*
*    Enable echoing of parameter file to log file as it is read. This
*    option is often useful for debugging the parameter file.
*
*  -r file
*
*   Set the program to reads a binary restart file to continue a
*   previous run.  The restart file name is given by the argument.
*
*  -p file
*
*   Set the parameter file name, given by the argument "file". 
*   The -p and -r options are incompatible. 
*
*  -c file
*
*   Set the command file name, given by the argument "file". The
*   command file may also be specified in the parameter file.
*
*  -i prefix
*
*   Set the input file path prefix, given by the argument "prefix".
*
*  -o prefix
*
*   Set the output file path prefix, given by the argument "prefix".
*
*  -f
*
*   Set replicated mode for parallel simulations.
*
* Input and output files:
*
* Serial: If compiled in serial mode, with MPI disabled (ifndef UTIL_MPI), 
* the default file Log::file() is written to standard output (std::cout).
* If the param file is not specified using a -p option, the param file 
* stream defaults to standard input. 
*
* Parallel: If compiled with MPI enabled (ifdef UTIL_MPI), paths for all 
* input output files associated with the processor with rank n are prefixed 
* by "n/". All output files produced by processor 6, for example, are thus
* put in a directory named "6/".  The default log file for processor n is
* is "n/log". 
*
* Parallel independent mode: If compiled with MPI enabled (ifdef UTIL_MPI) 
* and invoked without the -f option, each processor reads from a different 
* parameter file and command file which are in the numbered system directory 
* associated with that processor (or system).  If no param file name is 
* specified by a -p option, the param file name for processor n defaults 
* to "n/param". 
*
* Parallel replicated mode: If compiled with MPI enabled (ifdef UTIL_MPI) 
* and invoked with the "-f" option, a single parameter file and control 
* file is read from std::cin, as in a serial job, using a "perturbation"
* to set a series of slightly different parameters on different processors.
* If no param file name is specified by a -p option, the param file stream 
* defaults to standard input. 
*/

int main(int argc, char **argv)
{
   #ifdef UTIL_MPI
   MPI::Init();
   McMd::MdSimulation simulation(MPI::COMM_WORLD);
   #else
   McMd::MdSimulation simulation;
   #endif

   // Process command line options
   simulation.setOptions(argc, argv);

   // Read parameters from default parameter file
   simulation.readParam();

   // Read command script to run simulation
   simulation.readCommands();

   #ifdef UTIL_MPI
   if (MPI::Is_initialized()) {
      MPI::Finalize();
   }
   #endif

   return 0;
}
