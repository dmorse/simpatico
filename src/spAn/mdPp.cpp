/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <spAn/processor/Processor.h>

/**
* \page mdPp_page mdPp - postprocessing analysis program.
*
* Single-processor analysis program for postprocessing MD trajectories.
* 
* Usage:
*
*    mdPp [-e] param [command]
*
* Options:
*
*     -e enable echoing of the param file (no argument)
*
* Non-option arguments:
*
*     param   - name of parameter file
*
*     command - name of command file script (optional).
*
* If no command file is provided, commands are read from standard
* input (i.e., from the keyboard).
*
* \ingroup SpAn_Module
*/

int main(int argc, char** argv)
{
   SpAn::Processor processor;
   processor.setOptions(argc, argv);
   processor.readParam();
   processor.readCommands();
   return 0;
}
