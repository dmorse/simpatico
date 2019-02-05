#include "DomainTest.h"
int main()
{
   #ifdef UTIL_MPI 
   MPI_Init(&argc, &argv);
   IntVector::commitMpiType();
   #endif 

   TEST_RUNNER(DomainTest) runner1;
   runner1.run();

   #ifdef UTIL_MPI
   MPI_Finalize();
   #endif

} 
