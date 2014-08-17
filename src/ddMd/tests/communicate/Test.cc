#include "CommunicateTestComposite.h"

int main()
{
   #ifdef UTIL_MPI 
   MPI::Init();
   IntVector::commitMpiType();
   Vector::commitMpiType();

   CommunicateTestComposite runner;
   runner.run();

   MPI::Finalize();
   #endif

} 
