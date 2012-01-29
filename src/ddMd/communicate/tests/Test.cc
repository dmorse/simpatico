#include "DomainTest.h"
#include "BufferTest.h"
#include "AtomDistributorTest.h"
#include "ExchangerTest.h"

int main()
{
   #ifdef UTIL_MPI 
   MPI::Init();
   IntVector::commitMpiType();
   #endif 

   TEST_RUNNER(DomainTest) runner1;
   runner1.run();

   #ifdef UTIL_MPI
   TEST_RUNNER(BufferTest) runner2;
   runner2.run();

   TEST_RUNNER(AtomDistributorTest) runner3;
   runner3.run();

   TEST_RUNNER(ExchangerTest) runner4;
   runner4.run();

   MPI::Finalize();
   #endif

} 
