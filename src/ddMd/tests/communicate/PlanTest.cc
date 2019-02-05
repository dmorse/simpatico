#include "PlanTest.h"
int main()
{
   #ifdef UTIL_MPI 
   MPI_Init(&argc, &argv);
   #endif 

   TEST_RUNNER(PlanTest) runner1;
   runner1.run();

   #ifdef UTIL_MPI
   MPI_Finalize();
   #endif

} 
