/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <tools/processor/Processor.h>

/**
* \page mdPp_page mdPp - postprocessing analysis program
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
*
*      Enable echoing of the param file to std::out (no argument)
*
*     -p paramFileName
*
*      Specify name of a parameter file name as argument
*
*     -p commandFileName
*
*      Specify the name of a command file name as the argument
*
* If option -p is not set, so that no command file name is given,
* commands are read from standard input (i.e., from the keyboard).
*/

int main(int argc, char** argv)
{
   Tools::Processor processor;
   processor.setOptions(argc, argv);
   processor.readParam();
   processor.readCommands();
   return 0;
}
