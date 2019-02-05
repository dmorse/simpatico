#include "ExchangerForceTest.h"
int main()
{
   #ifdef UTIL_MPI
   MPI_Init(&argc, &argv);
   IntVector::commitMpiType();
   Vector::commitMpiType();

   TEST_RUNNER(ExchangerForceTest) runner;
   runner.run();

   MPI_Finalize();
   #endif

}

