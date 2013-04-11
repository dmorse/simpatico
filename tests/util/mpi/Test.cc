#include "MpiSendRecvTest.h"
#include "MpiFileIoTest.h"
#include "MpiLoggerTest.h"

using namespace Util;

int main()
{
   MPI::Init();
   Vector::commitMpiType();
   IntVector::commitMpiType();

   TEST_RUNNER(MpiSendRecvTest) test1;
   test1.run();

   TEST_RUNNER(MpiFileIoTest) test2;
   test2.run();

   //TEST_RUNNER(MpiLoggerTest) test3;
   //test3.run();

   MPI::Finalize();
}
