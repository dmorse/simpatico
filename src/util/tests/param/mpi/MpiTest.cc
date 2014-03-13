#ifdef UTIL_MPI

#define TEST_MPI
#include "MpiParamTestComposite.h"

#include <util/space/Vector.h>
#include <util/space/IntVector.h>
#include <util/global.h>

int main()
{
   MPI::Init();
   Vector::commitMpiType();
   IntVector::commitMpiType();

   MpiParamTestComposite runner;
   runner.run();

   MPI::Finalize();
}
#endif
