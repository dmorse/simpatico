//#include "CellTest.h"
//#include "CellListTest.h"
#include "NeighborTestComposite.h"

int main() 
{

   #ifdef UTIL_MPI 
   MPI_Init(&argc, &argv);
   //IntVector::commitMpiType();
   //Vector::commitMpiType();
   #endif

   NeighborTestComposite runner;
   runner.run();

   #ifdef UTIL_MPI 
   MPI_Finalize();
   #endif

}
