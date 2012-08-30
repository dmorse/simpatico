#include "GroupTest.h"
int main()
{
   #ifdef UTIL_MPI 
   MPI::Init();
   Vector::commitMpiType();
   #endif

   TEST_RUNNER(GroupTest) runner;
   runner.run();

   #ifdef UTIL_MPI 
   MPI::Finalize();
   #endif

} 
