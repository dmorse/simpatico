#include "AngleTest.h"
#include "CosineAngleTest.h"
#include "CosineSqAngleTest.h"
#include "HarmonicAngleTest.h"

int main() 
{
   TEST_RUNNER(AngleTest) runner1;
   runner1.run();

   TEST_RUNNER(CosineAngleTest) runner2;
   runner2.run();

   TEST_RUNNER(CosineSqAngleTest) runner3;
   runner3.run();

   TEST_RUNNER(HarmonicAngleTest) runner4;
   runner4.run();

}
