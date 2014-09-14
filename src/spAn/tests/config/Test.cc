#include "DdMdConfigReaderTest.h"
#include "HoomdConfigReaderTest.h"
#include "TypeMapTest.h"

int main()
{
   //StorageTestComposite runner;
   //runner.run();

   TEST_RUNNER(DdMdConfigReaderTest) runner1;
   runner1.run();

   TEST_RUNNER(TypeMapTest) runner2;
   runner2.run();

   TEST_RUNNER(HoomdConfigReaderTest) runner3;
   runner3.run();

} 
