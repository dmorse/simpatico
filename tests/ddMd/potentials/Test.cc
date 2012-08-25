#include "PairPotentialTest.h"

int main()
{

   #ifdef UTIL_MPI 
   MPI::Init();
   IntVector::commitMpiType();
   Vector::commitMpiType();
   #endif

   TEST_RUNNER(PairPotentialTest) runner;
   runner.run();

   #ifdef UTIL_MPI 
   MPI::Finalize();
   #endif

}

