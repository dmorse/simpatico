#include "PlanTest.h"
int main()
{
   #ifdef UTIL_MPI 
   MPI::Init();
   #endif 

   TEST_RUNNER(PlanTest) runner1;
   runner1.run();

   #ifdef UTIL_MPI
   MPI::Finalize();
   #endif

} 
