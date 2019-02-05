#include "ExchangerTest.h"
int main()
{
   #ifdef UTIL_MPI
   MPI_Init(&argc, &argv);
   IntVector::commitMpiType();
   Vector::commitMpiType();

   TEST_RUNNER(ExchangerTest) runner;
   runner.run();

   MPI_Finalize();
   #endif

}

