#include "ModifierTest.h"
#include "ModifierManagerTest.h"
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

   TEST_RUNNER(ModifierManagerTest) runner2;
   runner2.run();

   #ifdef UTIL_MPI 
   MPI::Finalize();
   #endif

}
