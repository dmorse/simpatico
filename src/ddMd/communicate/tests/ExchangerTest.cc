#include "ExchangerTest.h"
int main()
{
   #ifdef UTIL_MPI
   MPI::Init();
   IntVector::commitMpiType();

   TEST_RUNNER(ExchangerTest) runner;
   runner.run();

   MPI::Finalize();
   #endif

}

