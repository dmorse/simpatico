
#include <util/global.h>
#define TEST_MPI
#include "MpiParamTest.h"
#include "MpiParamCompositeTest.h"
#include "MpiManagerTest.h"
#include "MpiParamTestComposite.h"

#include <util/space/Vector.h>
#include <util/space/IntVector.h>

int main()
{
   MPI::Init();
   Vector::commitMpiType();
   IntVector::commitMpiType();

   //TEST_RUNNER(MpiParamTest) test1;
   //test1.run();

   //TEST_RUNNER(MpiParamCompositeTest) test2;
   //test2.run();

   //TEST_RUNNER(MpiManagerTest) test3;
   //test3.run();

   MpiParamTestComposite runner;
   runner.run();

   MPI::Finalize();
}
#ifdef UTIL_MPI
#endif
