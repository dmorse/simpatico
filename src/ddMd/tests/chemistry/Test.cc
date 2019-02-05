#include "ChemistryTestComposite.h"

int main()
{
   #ifdef UTIL_MPI 
   MPI_Init(&argc, &argv);
   #endif

   ChemistryTestComposite runner;
   runner.run();

   #ifdef UTIL_MPI
   MPI_Finalize();
   #endif

} 
