#include "AtomStorageTest.h"
#include "SystemTest.h"

int main()
{
   //SystemTestComposite runner;
   //runner.run();

   TEST_RUNNER(AtomStorageTest) runner1;
   runner1.run();

   TEST_RUNNER(SystemTest) runner2;
   runner2.run();
} 
