#include "BufferTest.h"

int main()
{
   #ifdef UTIL_MPI 
   MPI::Init();
   #endif

   TEST_RUNNER(BufferTest) runner;
   runner.run();

   #ifdef UTIL_MPI 
   MPI::Finalize();
   #endif

} 
