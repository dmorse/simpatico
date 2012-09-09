/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/mcSimulation/McSimulation.h>
#include "SliplinkMcModule.h"

/**
* Program for Monte Carlo simulation.
*
* Command line options:
*
*   -e  
*    Enable echoing of the parameter file as it is read.
*
*   -p  
*    Enable use of a free energy perturbation. 
*
*   -r filename
*    Restart a simulation. The parameter "filename" is the 
*    base name for the 3 restart files, named filename.prm, 
*    filename.rst, and filename.cmd.
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
   #ifdef UTIL_MPI
   MPI::Init();
   McMd::McSimulation simulation(MPI::COMM_WORLD);
   #else
   McMd::McSimulation simulation;
   #endif

   // Add Module
   McMd::SliplinkMcModule sliplinkMcModule(simulation);
   sliplinkMcModule.addFactories();

   // Process command line options
   simulation.setOptions(argc, argv);

   // Read parameters from default parameter file
   simulation.readParam();
   
   // Read command script to run simulation
   simulation.readCommands();

   #ifdef UTIL_MPI
   MPI::Finalize();
   #endif

   return 0;
}
