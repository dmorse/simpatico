#include "AtomCollectorTest.h"

int main()
{
   #ifdef UTIL_MPI
   MPI_Init(&argc, &argv);
   IntVector::commitMpiType();
   #endif

   TEST_RUNNER(AtomCollectorTest) runner;
   runner.run();

   #ifdef UTIL_MPI
   MPI_Finalize();
   #endif

}

