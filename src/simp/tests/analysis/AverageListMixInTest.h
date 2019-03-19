#ifndef SIMP_AVERAGE_LIST_MIXIN_TEST_H
#define SIMP_AVERAGE_LIST_MIXIN_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include "MockAverageList.h"
#include <simp/analysis/AverageListMixIn.h>
#include <util/accumulators/Average.h>
#include <util/archives/BinaryFileOArchive.h>
#include <util/misc/FileMaster.h>

#include <iostream>
#include <fstream>

using namespace Util;
using namespace Simp;

class AverageListMixInTest : public UnitTest 
{

private:

   FileMaster fileMaster;
   MockAverageList analyzer;

public:

   AverageListMixInTest() :
     fileMaster(),
     analyzer(fileMaster)
   {}

   void setUp() 
   { 
      //setVerbose(2);
      fileMaster.setOutputPrefix(filePrefix());

   } 

   void tearDown() 
   {}

   void testReadParam();

   void testSample();

   void testSaveLoad();
};

void AverageListMixInTest::testReadParam()
{
   printMethod(TEST_FUNC);

   std::ifstream in;
   openInputFile("in/MockAverageList", in);
   analyzer.readParam(in);
   in.close();

   if (verbose() > 1) {
      std::cout << std::endl;
      analyzer.writeParam(std::cout);
      std::cout << std::endl;
   }

}

void AverageListMixInTest::testSample()
{
   printMethod(TEST_FUNC);

   std::ifstream in;
   openInputFile("in/MockAverageList", in);
   analyzer.readParam(in);
   in.close();

   analyzer.setup();
   int i, j;
   j = 0;
   for (i = 0; i < 1000; ++i) {
      analyzer.sample(j);
      j += analyzer.interval();
   }
   analyzer.output();
}

void AverageListMixInTest::testSaveLoad()
{
   printMethod(TEST_FUNC);

   std::ifstream in;
   openInputFile("in/MockAverageList", in);
   analyzer.readParam(in);
   in.close();

   analyzer.setup();
   int i, j;
   j = 0;
   for (i = 0; i < 100; ++i) {
      analyzer.sample(j);
      j += analyzer.interval();
   }

   std::ofstream ofile;
   openOutputFile("tmp/arAveList", ofile);
   BinaryFileOArchive oar(ofile);
   analyzer.save(oar);
   ofile.close();

   std::ifstream ifile;
   openInputFile("tmp/arAveList", ifile);
   BinaryFileIArchive iar(ifile);
   MockAverageList analyzer2(fileMaster);
   analyzer2.load(iar);
   ifile.close();
    
   TEST_ASSERT(analyzer.nSamplePerBlock() == analyzer2.nSamplePerBlock());
   double ave0 = analyzer2.accumulator(0).average();
   double ave1 = analyzer2.accumulator(1).average();
   TEST_ASSERT(fabs(ave0) < 0.1);
   TEST_ASSERT(fabs(ave1 - 1.0) < 0.1);

   if (verbose() > 1) {
      std::cout << "\n";
      std::cout << "average(0) = " << ave0 << "\n";
      std::cout << "average(1) = " << ave1 << "\n";
   }

}

TEST_BEGIN(AverageListMixInTest)
TEST_ADD(AverageListMixInTest, testReadParam)
TEST_ADD(AverageListMixInTest, testSample)
#ifndef UTIL_MPI
TEST_ADD(AverageListMixInTest, testSaveLoad)
#endif
TEST_END(AverageListMixInTest)

#endif
