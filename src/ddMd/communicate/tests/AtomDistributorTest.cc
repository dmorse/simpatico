#include "AtomDistributorTest.h"

int main()
{
   #ifdef UTIL_MPI
   MPI::Init();
   IntVector::commitMpiType();
   #endif

   TEST_RUNNER(AtomDistributorTest) runner;
   runner.run();

   #ifdef UTIL_MPI
   MPI::Finalize();
   #endif

}

