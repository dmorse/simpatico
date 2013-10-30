#include "ModifierTest.h"
//#include "ModifierTestComposite.h"

int main() 
{

   #ifdef UTIL_MPI 
   MPI::Init();
   IntVector::commitMpiType();
   Vector::commitMpiType();
   #endif

   //ModifierTestComposite runner;
   TEST_RUNNER(ModifierTest) runner;
   runner.run();

   #ifdef UTIL_MPI 
   MPI::Finalize();
   #endif

}
