#include "DomainTest.h"
#include "BufferTest.h"
#include "AtomDistributorTest.h"
#include "ExchangerTest.h"

int main()
{
   #ifdef UTIL_MPI 
   MPI::Init();
   IntVector::commitMpiType();
   Vector::commitMpiType();

   TEST_RUNNER(DomainTest) runner1;
   runner1.run();

   TEST_RUNNER(BufferTest) runner2;
   runner2.run();

   TEST_RUNNER(AtomDistributorTest) runner3;
   runner3.run();

   TEST_RUNNER(AtomDistributorTest) runner4;
   runner4.run();

   TEST_RUNNER(ExchangerTest) runner5;
   runner5.run();

   MPI::Finalize();
   #endif

} 
