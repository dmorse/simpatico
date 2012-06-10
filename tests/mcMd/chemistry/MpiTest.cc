#ifdef UTIL_MPI
#include "MpiChemistryTest.h"

using namespace Util;

int main()
{
   MPI::Init();
   McMd::commitMpiTypes();

   TEST_RUNNER(MpiChemistryTest) test;
   test.run();

   MPI::Finalize();
}
#endif
