#include "SerializeConfigIoTest.h"

int main()
{
   #ifdef UTIL_MPI
   MPI::Init();
   IntVector::commitMpiType();
   Vector::commitMpiType();
   #endif

   TEST_RUNNER(SerializeConfigIoTest) runner;
   runner.run();

   #ifdef UTIL_MPI
   MPI::Finalize();
   #endif

}

