#include "CommunicateTestComposite.h"

int main(int argc, char** argv)
{
   #ifdef UTIL_MPI 
   MPI_Init(&argc, &argv);
   IntVector::commitMpiType();
   Vector::commitMpiType();

   CommunicateTestComposite runner;
   runner.run();

   MPI_Finalize();
   #endif

} 
