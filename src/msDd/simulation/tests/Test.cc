
//#include <util/space/Vector.h>
//#include <util/space/IntVector.h>

#include <util/util/initStatic.h>
#include <util/global.h>
#define TEST_MPI
#include "SimulationTest.h"

int main()
{
   #ifdef UTIL_MPI 
   MPI::Init();
   Util::initStatic();
   //Util::IntVector::commitMpiType();
   //Util::Vector::commitMpiType();
   #endif 

   TEST_RUNNER(SimulationTest) runner;
   runner.run();

   #ifdef UTIL_MPI
   if (MPI::Is_initialized()) {
      MPI::Finalize();
   }
   #endif

} 
