#include "SpeciesTestComposite.h"
//#include "SpeciesTest.h"
//#include "HomopolymerTest.h"
//#include "MultiblockTest.h"
//#include "HomoRingTest.h"

int main()
{
   SpeciesTestComposite test;
   test.run();

   #if 0
   TEST_RUNNER(SpeciesTest) test1;
   test1.run();

   TEST_RUNNER(HomopolymerTest) test2;
   test2.run();

   TEST_RUNNER(HomoRingTest) test3;
   test3.run();

   TEST_RUNNER(MultiblockTest) test4;
   test4.run();
   #endif
}
