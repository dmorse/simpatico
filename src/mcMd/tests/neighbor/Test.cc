
#ifdef UTIL_MPI
#define TEST_MPI
#endif

#include "NeighborTestComposite.h"

int main(int argc, char** argv) 
{
   #ifdef UTIL_MPI
   MPI_Init(&argc, &argv);
   #endif

   NeighborTestComposite runner;
   runner.run();

   #ifdef UTIL_MPI
   MPI_Finalize();
   #endif
}
