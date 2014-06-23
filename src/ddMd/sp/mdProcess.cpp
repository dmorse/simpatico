/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/sp/processor/Processor.h>

// std headers
#include <fstream>
#include <unistd.h>
#include <stdlib.h>

using namespace Util;

   /**
   * Program for postprocessing ddSim MD trajectories.
   */
   int main(int argc, char** argv)
   {

      DdMd::Processor processor;

      // Read command-line arguments
      bool tFlag = false;
      int c;
      opterr = 0;
      while ((c = getopt(argc, argv, "t")) != -1) {
         switch (c) {
         case 't':
           tFlag = true;
           break;
         case '?':
           Log::file() << "Unknown option -" << optopt << std::endl;
         }
      }
      if (argc - optind != 4) {
        std::cout << "optind = " << optind << std::endl;
        std::cout << "argc   = " << argc << std::endl;
         UTIL_THROW("Wrong number of arguments");
      }
      const char* paramFileName = argv[optind];
      const char* configFileName = argv[optind+1];
      int first = atoi(argv[optind+2]);
      int last = atoi(argv[optind+3]);
      std::cout << "paramFile  = " << paramFileName << std::endl;
      std::cout << "configFile = " << configFileName << std::endl;
      std::cout << "first      = " << first << std::endl;
      std::cout << "last       = " << last << std::endl;
   
      // Read parameter file 
      std::ifstream paramFile;
      paramFile.open(paramFileName); 
      processor.readParam(paramFile);
      paramFile.close();

      // Process dumps
      processor.analyzeDumps(first, last, paramFileName);

      return 0;
   }

