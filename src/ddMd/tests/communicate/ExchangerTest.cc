#include "ExchangerTest.h"
int main(int argc, char** argv)
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

