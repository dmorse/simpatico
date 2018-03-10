#include "ToolsTestComposite.h"

/*
* This program runs all unit tests in the tools directory.
*/ 
int main(int argc, char* argv[])
{
   try {

      ToolsTestComposite runner;
      if (argc > 2) {
         UTIL_THROW("Too many arguments");
      }
      if (argc == 2) {
         runner.addFilePrefix(argv[1]);
      }

      // Run unit tests, count failures
      int failures = runner.run();

      return (failures != 0);
   } catch (...) {
      std::cerr << "Uncaught exception in tools Test.cc" 
                << std::endl;
      return 1;
   }

}
