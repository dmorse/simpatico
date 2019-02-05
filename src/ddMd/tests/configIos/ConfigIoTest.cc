#include "ConfigIoTest.h"

int main()
{
   #ifdef UTIL_MPI
   MPI_Init(&argc, &argv);
   IntVector::commitMpiType();
   Vector::commitMpiType();
   #endif

   TEST_RUNNER(ConfigIoTest) runner;
   runner.run();

   #ifdef UTIL_MPI
   MPI_Finalize();
   #endif

}

