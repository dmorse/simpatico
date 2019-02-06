#include "SimulationTest.h"
//#include "SimulationTest2.h"
int main(int argc, char** argv)
{
   #ifdef UTIL_MPI 
   MPI_Init(&argc, &argv);
   IntVector::commitMpiType();
   Vector::commitMpiType();
   #endif 

   TEST_RUNNER(SimulationTest) runner;
   runner.run();

   #ifdef UTIL_MPI
   MPI_Finalize();
   #endif

} 
