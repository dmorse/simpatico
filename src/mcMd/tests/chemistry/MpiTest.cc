#ifdef UTIL_MPI
#define TEST_MPI
#include "MpiChemistryTest.h"

using namespace Util;

int main(int argc, char** argv)
{
   MPI_Init(&argc, &argv);
   McMd::commitMpiTypes();

   TEST_RUNNER(MpiChemistryTest) test;
   test.run();

   MPI_Finalize();
}
#endif
