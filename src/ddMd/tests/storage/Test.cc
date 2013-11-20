//#include "AtomMapTest.h"
//#include "AtomStorageTest.h"
//#include "BondStorageTest.h"
#include "StorageTestComposite.h"

int main()
{
   #ifdef UTIL_MPI 
   MPI::Init();
   #endif

   StorageTestComposite runner;
   runner.run();

   #if 0
   TEST_RUNNER(AtomMapTest) runner1;
   runner1.run();

   TEST_RUNNER(AtomStorageTest) runner2;
   runner2.run();

   TEST_RUNNER(BondStorageTest) runner3;
   runner3.run();
   #endif

   #ifdef UTIL_MPI
   MPI::Finalize();
   #endif

} 
