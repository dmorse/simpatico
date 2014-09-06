/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <spAn/processor/Processor.h>

// std headers
#include <fstream>
#include <unistd.h>
#include <stdlib.h>

using namespace Util;

/**
* \page mdPp_page mdPp - postprocessing analysis program.
*
* Single-processor analysis program for postprocessing MD trajectories.
* 
* Usage:
*
*    mdPp [-t] [-i inputConfigReader]  param config first last
*
* Required arguments:
*
*     param  - name of parameter file
*
*     config - base name of configuration file(s)
*
*     first  - index of first configuration or frame
*
*     last   - index of last configuration or frame
*
* Options:
*
*     -t
*
*        Sets to read a trajectory file. If not set, reads a sequence
*        of configuration files.
*
*     -i inputConfigReader
*
*        Specify input configuration/trajectory format. The required
*        argument inputConfigReader is the name of a ConfigReader subclass.
*
* \ingroup SpAn_Module
*/

int main(int argc, char** argv)
{

   SpAn::Processor processor;

   // Parse command-line arguments
   bool tFlag = false;
   bool cFlag = false;
   bool dFlag = false;
   std::string configStyle;
   std::string trajectoryStyle;
   int c;
   opterr = 0;
   while ((c = getopt(argc, argv, "tc:d:")) != -1) {
      switch (c) {
      case 't':
         tFlag = true;
         break;
      case 'c':
         cFlag = true;
         configStyle = optarg;
         break;
      case 'd':
         dFlag = true;
         trajectoryStyle = optarg;
         break;
      case '?':
         Log::file() << "Unknown option -" << optopt << std::endl;
      }
   }

   // Read required arguments
   if (tFlag) {
      if (argc - optind != 3) {
         Log::file() << "optind = " << optind << std::endl;
         Log::file() << "argc   = " << argc << std::endl;
         UTIL_THROW("Wrong number of required arguments");
      }
   } else {
      if (argc - optind != 4) {
         Log::file() << "optind = " << optind << std::endl;
         Log::file() << "argc   = " << argc << std::endl;
         UTIL_THROW("Wrong number of required arguments");
      }
   }

   // Read parameter file
   const char* paramFileName = argv[optind];
   Log::file() << "paramFile  = " << paramFileName << std::endl;
   std::ifstream paramFile;
   paramFile.open(paramFileName);
   processor.readParam(paramFile);
   paramFile.close();

   const char* configFileName = argv[optind+1];
   Log::file() << "configFile = " << configFileName << std::endl;

   if (cFlag) {
      Log::file() << "configStyle= " 
                  << configStyle << std::endl;
      processor.setConfigReader(configStyle);
   }

   if (tFlag) {
       processor.readConfig(configFileName);
       const char* trajectoryFileName = argv[optind+2];
       Log::file() << "trajectoryFile = " << trajectoryFileName << std::endl;
       if (dFlag) {
          Log::file() << "trajectoryStyle= " 
                      << trajectoryStyle << std::endl;
          processor.setTrajectoryReader(trajectoryStyle);
       }
       processor.analyzeTrajectory(trajectoryFileName);
   } else {
      int first = atoi(argv[optind+2]);
      int last = atoi(argv[optind+3]);
      Log::file() << "first = " << first << std::endl;
      Log::file() << "last  = " << last << std::endl;
      processor.analyzeDumps(first, last, configFileName);
   }



   // Process dumps

   return 0;
}

