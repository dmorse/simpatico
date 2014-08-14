#ifndef SPAN_TEST_CPP
#define SPAN_TEST_CPP
#include "SpAnTestComposite.h"

/*
* This program runs all unit tests in the spAn directory.
*/ 
int main(int argc, char* argv[])
{
   SpAnTestComposite runner;
   if (argc > 2) {
      UTIL_THROW("Too many arguments");
   }
   if (argc == 2) {
      runner.addFilePrefix(argv[1]);
   }
   runner.run();
}
#endif 
