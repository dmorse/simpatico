#include "DdMdConfigReaderTest.h"
#include "HoomdConfigReaderTest.h"

int main()
{
   //StorageTestComposite runner;
   //runner.run();

   TEST_RUNNER(DdMdConfigReaderTest) runner1;
   runner1.run();

   TEST_RUNNER(HoomdConfigReaderTest) runner2;
   runner2.run();
} 
