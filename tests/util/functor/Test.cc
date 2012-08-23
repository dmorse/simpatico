#include "MethodFunctorTest.h"

int main() 
{

   #ifdef UTIL_MPI
   MPI::Init();
   #endif

   TEST_RUNNER(MethodFunctorTest) test;
   test.run();

   #ifdef UTIL_MPI
   MPI::Finalize();
   #endif

}
