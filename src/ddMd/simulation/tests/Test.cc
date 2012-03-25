#include "SystemTest.h"
//#include "SystemTest2.h"
int main()
{
   #ifdef UTIL_MPI 
   MPI::Init();
   IntVector::commitMpiType();
   Vector::commitMpiType();
   #endif 

   TEST_RUNNER(SystemTest) runner;
   runner.run();

   #ifdef UTIL_MPI
   MPI::Finalize();
   #endif

} 
