#include "AtomTest.h"
int main()
{
   #ifdef UTIL_MPI 
   MPI_Init(&argc, &argv);
   Vector::commitMpiType();
   #endif

   TEST_RUNNER(AtomTest) runner;
   runner.run();

   #ifdef UTIL_MPI 
   MPI_Finalize();
   #endif

} 
