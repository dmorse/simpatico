#ifndef SIMP_AVERAGE_MIXIN_TEST_H
#define SIMP_AVERAGE_MIXIN_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include "MockAverage.h"
#include <simp/analysis/AverageMixIn.h>
#include <util/misc/FileMaster.h>
#include <util/accumulators/Average.h>
#include <util/archives/BinaryFileOArchive.h>
#include <util/archives/BinaryFileIArchive.h>

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
   } 

   void tearDown() 
   {}

   void testReadParam();

   void testSample();

   void testSaveLoad();

};

void AverageMixInTest::testReadParam()
{
   printMethod(TEST_FUNC);

   std::ifstream in;
   openInputFile("in/MockAverage", in);
   analyzer.readParam(in);
   in.close();

   if (verbose() > 1) {
      std::cout << std::endl;
      analyzer.writeParam(std::cout);
      std::cout << std::endl;
   }

}

void AverageMixInTest::testSample()
{
   printMethod(TEST_FUNC);

   std::ifstream in;
   openInputFile("in/MockAverage", in);
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

void AverageMixInTest::testSaveLoad()
{
   printMethod(TEST_FUNC);

   std::ifstream in;
   openInputFile("in/MockAverage", in);
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
   openOutputFile("tmp/arAve", ofile);
   BinaryFileOArchive oar(ofile);
   analyzer.save(oar);
   ofile.close();

   std::ifstream ifile;
   openInputFile("tmp/arAve", ifile);
   BinaryFileIArchive iar(ifile);

   MockAverage analyzer2(fileMaster);
   analyzer2.load(iar);
   ifile.close();
    
   TEST_ASSERT(analyzer.nSamplePerBlock() == analyzer2.nSamplePerBlock());
   double ave = analyzer2.accumulator().average();
   TEST_ASSERT(fabs(ave - 1.0) < 0.1);

   if (verbose() > 1) {
      std::cout << "\n";
      std::cout << ave << "\n";
   }

}


TEST_BEGIN(AverageMixInTest)
TEST_ADD(AverageMixInTest, testReadParam)
TEST_ADD(AverageMixInTest, testSample)
#ifndef UTIL_MPI
TEST_ADD(AverageMixInTest, testSaveLoad)
#endif
TEST_END(AverageMixInTest)

#endif
