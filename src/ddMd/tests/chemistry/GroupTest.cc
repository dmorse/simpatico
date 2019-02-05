#include "GroupTest.h"
int main()
{
   #ifdef UTIL_MPI 
   MPI_Init(&argc, &argv);
   Vector::commitMpiType();
   #endif

   TEST_RUNNER(GroupTest) runner;
   runner.run();

   #ifdef UTIL_MPI 
   MPI_Finalize();
   #endif

} 
