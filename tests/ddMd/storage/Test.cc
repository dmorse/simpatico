#include "AtomStorageTest.h"
#include "BondStorageTest.h"

int main()
{
   #ifdef UTIL_MPI 
   MPI::Init();
   #endif

   TEST_RUNNER(AtomStorageTest) runner1;
   runner1.run();

   TEST_RUNNER(BondStorageTest) runner2;
   runner2.run();

   #ifdef UTIL_MPI
   MPI::Finalize();
   #endif

} 
