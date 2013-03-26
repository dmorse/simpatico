#include "AtomTest.h"
int main()
{
   #ifdef UTIL_MPI 
   MPI::Init();
   Vector::commitMpiType();
   #endif

   TEST_RUNNER(AtomTest) runner;
   runner.run();

   #ifdef UTIL_MPI 
   MPI::Finalize();
   #endif

} 
