#include "AtomStorageTest.h"
#include "StorageTest.h"

int main()
{
   //StorageTestComposite runner;
   //runner.run();

   TEST_RUNNER(AtomStorageTest) runner1;
   runner1.run();

   TEST_RUNNER(StorageTest) runner2;
   runner2.run();
} 
