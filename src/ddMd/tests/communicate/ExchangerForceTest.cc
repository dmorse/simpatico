#include "ExchangerForceTest.h"
int main()
{
   #ifdef UTIL_MPI
   MPI::Init();
   IntVector::commitMpiType();
   Vector::commitMpiType();

   TEST_RUNNER(ExchangerForceTest) runner;
   runner.run();

   MPI::Finalize();
   #endif

}

