#ifndef SIMP_AVERAGE_MIXIN_TEST_H
#define SIMP_AVERAGE_MIXIN_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include "MockAverage.h"
#include <simp/analysis/AverageMixIn.h>
#include <util/misc/FileMaster.h>

#include <iostream>
#include <fstream>

using namespace Util;
using namespace Simp;

class AverageMixInTest : public UnitTest 
{

private:

   FileMaster fileMaster;
   MockAverage analyzer;

public:

   AverageMixInTest() :
     fileMaster(),
     analyzer(fileMaster)
   {}

   void setUp() 
   { 
      //setVerbose(2);
      fileMaster.setOutputPrefix(filePrefix());

      std::ifstream in;
      openInputFile("in/MockAverage", in);
      analyzer.readParam(in);
      in.close();
   } 

   void tearDown() 
   {}

   void testReadParam();

   void testSample();

};

void AverageMixInTest::testReadParam()
{
   printMethod(TEST_FUNC);

   if (verbose() > 1) {
      std::cout << std::endl;
      analyzer.writeParam(std::cout);
      std::cout << std::endl;
   }

}

void AverageMixInTest::testSample()
{
   printMethod(TEST_FUNC);
   analyzer.setup();
   int i, j;
   j = 0;
   for (i = 0; i < 1000; ++i) {
      analyzer.sample(j);
      j += analyzer.interval();
   }
   analyzer.output();
}

TEST_BEGIN(AverageMixInTest)
TEST_ADD(AverageMixInTest, testReadParam)
TEST_ADD(AverageMixInTest, testSample)
TEST_END(AverageMixInTest)

#endif
