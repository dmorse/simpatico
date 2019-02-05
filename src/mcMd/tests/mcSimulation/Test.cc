#ifdef UTIL_MPI
#define TEST_MPI
#endif

#include "McSimulationTest.h"

int main(int argc, char** argv)
{
   #ifdef UTIL_MPI
   MPI_Init(&argc, &argv);
   #endif 

   TEST_RUNNER(McSimulationTest) runner;
   runner.run();

   #ifdef UTIL_MPI
   MPI_Finalize();
   #endif 
} 
