/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/mcSimulation/McSimulation.h>

/**
* \page mcSim_page mcSim - serial Monte Carlo program
*
* Single-processor Monte Carlo simulation program.
*
* Usage (single processor version, without MPI):
*
*    mcSim [-e] [-p] [-r restartFile] < paramFile
*
* Options:
*
*   -e  
*
*    Enable echoing of parameter file to log file as it is read. This
*    option is often useful for debugging the parameter file.
*
*   -p  
*
*    Enable use of a free energy perturbation. 
*
*  -r restartFile
*
*   This option reads a restart file and restarts a previous run. The 
*   command line argument restarFile argument is the shared base name of
*   the restart file and a corresponding command script file. The name
*   of the restart file is obtained by appending the file extension .rst 
*   to the base name, giving a name of the form restartFile.rst. The
*   corresponding command file must have the same base name and a file
*   name extension .cmd, giving a name of the form restartFile.cmd.
*
* Input and output files:
*
* If compiled in serial mode, with MPI disabled (ifndef UTIL_MPI), 
* the default parameter stream read by readParam() is read from 
* standard input (std::cin) and the default log file Log::file() 
* is written to standard output (std::cout).
*
* If compiled with MPI enabled (ifdef UTIL_MPI) and invoked with the 
* "-p" option, a single parameter file is read from std::cin, as in 
* a serial job.  
*
* If compiled with MPI enabled (ifdef UTIL_MPI) and invoked without 
* the -p option, each processor reads from a different parameter file, 
* which is named "n/param" for processor n.
*
* If compiled with MPI enabled (ifdef UTIL_MPI), paths for all output
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
