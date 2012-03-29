#include "GroupDistributorTest.h"

int main()
{
  #ifdef UTIL_MPI
   MPI::Init();
   IntVector::commitMpiType();
   Vector::commitMpiType();
   #endif

   TEST_RUNNER(GroupDistributorTest) runner;
   runner.run();

   #ifdef UTIL_MPI
   MPI::Finalize();
   #endif

}
