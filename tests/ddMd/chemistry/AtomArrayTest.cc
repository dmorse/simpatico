#include "AtomArrayTest.h"
int main()
{
   #ifdef UTIL_MPI 
   MPI::Init();
   Vector::commitMpiType();
   #endif

   TEST_RUNNER(AtomArrayTest) runner;
   runner.run();

   #ifdef UTIL_MPI 
   MPI::Finalize();
   #endif

} 
