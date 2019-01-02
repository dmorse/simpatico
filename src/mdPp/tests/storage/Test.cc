#include "AtomStorageTest.h"
#include "ConfigurationTest.h"

int main()
{
   //ConfigurationTestComposite runner;
   //runner.run();

   TEST_RUNNER(AtomStorageTest) runner1;
   runner1.run();

   TEST_RUNNER(ConfigurationTest) runner2;
   runner2.run();
} 
