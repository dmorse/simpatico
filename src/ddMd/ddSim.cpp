/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/simulation/Simulation.h>
#include <util/param/ParamComponent.h>
#include <util/global.h>

#include <iostream>

using namespace Util;

/*
* Main program for parallel domain-decomposition program (ddSim).
*/
int main(int argc, char **argv)
{

   #ifdef UTIL_MPI
   MPI_Init(&argc, &argv);
   #endif

   #ifdef UTIL_MPI
   DdMd::Simulation simulation(MPI_COMM_WORLD);
   #else
   DdMd::Simulation simulation();
   #endif

   // Read and apply command line options.
   simulation.setOptions(argc, argv); 

   // Read parameter file from default param stream (read on master).
   simulation.readParam(); 

   // Read command file (command file is specified in parameter file).
   simulation.readCommands();

   #ifdef UTIL_MPI
   MPI_Finalize();
   #endif

   return 0;
}
