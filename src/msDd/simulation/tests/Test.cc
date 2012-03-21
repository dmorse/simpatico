#include "SimulationTest.h"
//#include <util/space/Vector.h>
//#include <util/space/IntVector.h>

int main()
{
   #ifdef UTIL_MPI 
   MPI::Init();
   //Util::IntVector::commitMpiType();
   //Util::Vector::commitMpiType();
   #endif 

   TEST_RUNNER(SimulationTest) runner;
   runner.run();

   #ifdef UTIL_MPI
   MPI::Finalize();
   #endif

} 
