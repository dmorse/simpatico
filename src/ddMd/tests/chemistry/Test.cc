#include "ChemistryTestComposite.h"

int main()
{
   #ifdef UTIL_MPI 
   MPI::Init();
   #endif

   ChemistryTestComposite runner;
   runner.run();

   #ifdef UTIL_MPI
   MPI::Finalize();
   #endif

} 
