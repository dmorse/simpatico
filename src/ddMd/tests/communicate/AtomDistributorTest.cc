#include "AtomDistributorTest.h"

int main(int argc, char** argv)
{
   #ifdef UTIL_MPI
   MPI_Init(&argc, &argv);
   IntVector::commitMpiType();
   #endif

   TEST_RUNNER(AtomDistributorTest) runner;
   runner.run();

   #ifdef UTIL_MPI
   MPI_Finalize();
   #endif

}

