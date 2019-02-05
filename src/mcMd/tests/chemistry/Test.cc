#ifdef UTIL_MPI
#define TEST_MPI
#endif

#include "ChemistryTestComposite.h"
#include <util/misc/initStatic.h>

int main(int argc, char** argv) 
{
   #ifdef UTIL_MPI
   MPI_Init(&argc, &argv);
   #endif

   Util::initStatic();
   ChemistryTestComposite runner;
   runner.run();

   #ifdef UTIL_MPI
   MPI_Finalize();
   #endif
}
