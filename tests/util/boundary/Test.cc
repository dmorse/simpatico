//#include "BoundaryTestComposite.h"
#include "MonoclinicBoundaryTest.h"
#include "MonoclinicBoundaryMITest.h"

int main() 
{
   //BoundaryTestComposite runner;
   TEST_RUNNER(MonoclinicBoundaryMITest) runner;
   runner.run();

   return 0;
}
