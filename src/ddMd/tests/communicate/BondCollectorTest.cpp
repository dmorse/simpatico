#include "BondCollectorTest.h"

int main()
{
  #ifdef UTIL_MPI
   MPI::Init();
   IntVector::commitMpiType();
   Vector::commitMpiType();
   #endif

   TEST_RUNNER(BondCollectorTest) runner;
   runner.run();

   #ifdef UTIL_MPI
   MPI::Finalize();
   #endif

}
