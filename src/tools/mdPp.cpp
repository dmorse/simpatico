/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <tools/processor/Processor.h>

/**
* \page mdPp_page mdPp - postprocessing analysis program.
*
* Single-processor analysis program for postprocessing MD trajectories.
* 
* Usage:
*
*    mdPp [-e] -p param [-c command]
*
* Options:
*
*     -e 
*      Enable echoing of the param file (no argument)
*
*     -p paramFileName
*      Specify parameter file name as an argument
*
*     -p commandFileName
*      Specify command file name as an argument
*
* If option -p is not set, so that no command file name is given,
* commands are read from standard input (i.e., from the keyboard).
*
* \ingroup Tools_Module
*/

int main(int argc, char** argv)
{
   Tools::Processor processor;
   processor.setOptions(argc, argv);
   processor.readParam();
   processor.readCommands();
   return 0;
}
