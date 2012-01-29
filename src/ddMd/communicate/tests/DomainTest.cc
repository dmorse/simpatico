#include "DomainTest.h"
int main()
{
   #ifdef UTIL_MPI 
   MPI::Init();
   IntVector::commitMpiType();
   #endif 

   TEST_RUNNER(DomainTest) runner1;
   runner1.run();

   #ifdef UTIL_MPI
   MPI::Finalize();
   #endif

} 
