//#include "BoundaryTestComposite.h"
#include "MonoclinicBoundaryTest.h"

int main() 
{
   //BoundaryTestComposite runner;
   TEST_RUNNER(MonoclinicBoundaryTest) runner;
   runner.run();

   return 0;
}
