#ifndef TOOLS_TEST_CPP
#define TOOLS_TEST_CPP
#include "ToolsTestComposite.h"

/*
* This program runs all unit tests in the tools directory.
*/ 
int main(int argc, char* argv[])
{
   ToolsTestComposite runner;
   if (argc > 2) {
      UTIL_THROW("Too many arguments");
   }
   if (argc == 2) {
      runner.addFilePrefix(argv[1]);
   }
   runner.run();
}
#endif 
