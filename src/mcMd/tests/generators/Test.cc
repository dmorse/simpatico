#include "PointGeneratorTest.h"
#include "LinearGeneratorTest.h"

int main()
{ 
   TEST_RUNNER(PointGeneratorTest) runner1;
   runner1.run();

   TEST_RUNNER(LinearGeneratorTest) runner2;
   runner2.run();
} 
