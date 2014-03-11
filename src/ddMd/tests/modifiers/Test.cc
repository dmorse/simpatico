#include "ModifierTest.h"
// #include "ModifierManagerTest.h"
// #include "ModifierTestComposite.h"

using namespace Util;
using namespace DdMd;

int main() 
{

   #ifdef UTIL_MPI 
   #ifdef TEST_MPI
   MPI::Init();
   //IntVector::commitMpiType();
   //Vector::commitMpiType();
   #endif
   #endif

   //ModifierTestComposite runner;
   TEST_RUNNER(ModifierTest) runner;
   runner.run();

   // TEST_RUNNER(ModifierManagerTest) runner2;
   // runner2.run();

   #ifdef UTIL_MPI 
   #ifdef TEST_MPI
   MPI::Finalize();
   #endif
   #endif

}
