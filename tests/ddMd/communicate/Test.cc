#include "DomainTest.h"
#include "BufferTest.h"
#include "AtomDistributorTest.h"
#include "GroupDistributorTest.h"
#include "ExchangerTest.h"
#include "ExchangerForceTest.h"
#include "AtomCollectorTest.h"
#include "BondCollectorTest.h"

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

   TEST_RUNNER(GroupDistributorTest) runner4;
   runner4.run();

   TEST_RUNNER(ExchangerTest) runner5;
   runner5.run();

   //TEST_RUNNER(ExchangerForceTest) runner6;
   //runner6.run();

   TEST_RUNNER(AtomCollectorTest) runner7;
   runner7.run();

   TEST_RUNNER(BondCollectorTest) runner8;
   runner8.run();

   MPI::Finalize();
   #endif

} 
