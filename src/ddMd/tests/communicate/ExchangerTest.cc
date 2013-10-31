#include "ExchangerTest.h"
int main()
{
   #ifdef UTIL_MPI
   MPI::Init();
   IntVector::commitMpiType();
   Vector::commitMpiType();

   TEST_RUNNER(ExchangerTest) runner;
   runner.run();

   MPI::Finalize();
   #endif

}

