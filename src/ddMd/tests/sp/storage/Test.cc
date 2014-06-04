#include "SpAtomStorageTest.h"
#include "SpConfigurationTest.h"

int main()
{
   //SpConfigurationTestComposite runner;
   //runner.run();

   TEST_RUNNER(SpAtomStorageTest) runner1;
   runner1.run();

   TEST_RUNNER(SpConfigurationTest) runner2;
   runner2.run();
} 
