#ifdef UTIL_MPI
#define TEST_MPI
#endif

#include "SimulationTestComposite.h"

int main(int argc, char** argv)
{
   #if UTIL_MPI
   MPI_Init(&argc, &argv);
   #endif 

   SimulationTestComposite runner;
   runner.run();

   #if UTIL_MPI
   MPI_Finalize();
   #endif 
} 
