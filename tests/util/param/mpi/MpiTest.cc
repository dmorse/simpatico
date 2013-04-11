
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

   MpiParamTestComposite runner;
   runner.run();

   MPI::Finalize();
}
#ifdef UTIL_MPI
#endif
