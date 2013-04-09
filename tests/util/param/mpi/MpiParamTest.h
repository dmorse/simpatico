#ifndef MPI_PARAM_TEST_H
#define MPI_PARAM_TEST_H

#include <util/global.h>
#ifdef UTIL_MPI

#include <util/param/ScalarParam.h>

#ifndef TEST_MPI
#define TEST_MPI
#endif

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <iostream>
#include <fstream>

using namespace Util;

class MpiParamTest : public UnitTest
{

public:

   MpiParamTest()
    : UnitTest()
   {}

   void testScalarParamInt() 
   {
      printMethod(TEST_FUNC);

      std::ifstream in;
      int value;

      ScalarParam<int> param("MyLabel", value);
      param.setIoCommunicator(communicator());
      if (mpiRank() == 0)
         in.open("in/ScalarParamInt");
      param.readParam(in);
      if (mpiRank() == 1) {
         std::cout << std::endl;
         param.writeParam(std::cout);
         std::cout.flush();
      }
   }

};

TEST_BEGIN(MpiParamTest)
TEST_ADD(MpiParamTest, testScalarParamInt)
TEST_END(MpiParamTest)

#endif
#endif
