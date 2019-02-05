#include "BufferTest.h"

int main()
{
   #ifdef UTIL_MPI 
   MPI_Init(&argc, &argv);
   #endif

   TEST_RUNNER(BufferTest) runner;
   runner.run();

   #ifdef UTIL_MPI 
   MPI_Finalize();
   #endif

} 
