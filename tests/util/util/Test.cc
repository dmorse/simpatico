#include "ExceptionTest.h"
#include "ioUtilTest.h"

int main() 
{

   #ifdef UTIL_MPI
   MPI::Init();
   #endif

   TEST_RUNNER(ExceptionTest) test1;
   test1.run();

   TEST_RUNNER(ioUtilTest) test2;
   test2.run();

   #ifdef UTIL_MPI
   MPI::Finalize();
   #endif
}
