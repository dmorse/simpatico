#include "BufferTest.h"

int main(int argc, char** argv)
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
