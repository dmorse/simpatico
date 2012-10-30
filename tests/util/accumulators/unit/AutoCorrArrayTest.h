#ifndef AUTOCORR_ARRAY_TEST_H
#define AUTOCORR_ARRAY_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <util/accumulators/AutoCorrArray.h>
#include <util/archives/MemoryOArchive.h>
#include <util/archives/MemoryIArchive.h>
#include <util/archives/MemoryCounter.h>

#include <iostream>
#include <fstream>

using namespace Util;

class AutoCorrArrayTest : public UnitTest
{

public:

   AutoCorrArray<double, double> accumulator_;

   AutoCorrArrayTest() 
   {}

   void setUp() 
   { 
      std::ifstream paramFile;
      openInputFile("in/AutoCorrArray", paramFile); 
      accumulator_.readParam(paramFile);
      paramFile.close();
      accumulator_.setNEnsemble(4);
   }

   void tearDown() 
   {}

   void readData() 
   {
      DArray<double> data;
      int i, j, m, n, p;

      std::ifstream dataFile;
      openInputFile("in/data", dataFile); 
      dataFile >> m;
      n = accumulator_.nEnsemble();
      data.allocate(n);
      p = m / n;
      for (i = 0; i < p; ++i) {
         for (j = 0; j < n; ++j) {
           dataFile >> data[j];
         }
         accumulator_.sample(data);
      }
      dataFile.close();
   }

   void testReadParam() 
   {
      printMethod(TEST_FUNC);

      printEndl();
      accumulator_.writeParam(std::cout);
   }

   void testSample() 
   {
      printMethod(TEST_FUNC);

      readData();

      printEndl();
      accumulator_.output(std::cout);
   }

   void testSerialize() 
   {
      printMethod(TEST_FUNC);
      printEndl();

      readData();

      int size = memorySize(accumulator_);

      MemoryOArchive u;
      u.allocate(size);

      std::cout << size << std::endl;

      u << accumulator_;
      TEST_ASSERT(u.cursor() == u.begin() + size);

      MemoryIArchive v;
      v = u;

      AutoCorrArray<double, double> clone;
      v & clone;

      clone.output(std::cout);
   }

};

TEST_BEGIN(AutoCorrArrayTest)
TEST_ADD(AutoCorrArrayTest, testReadParam)
TEST_ADD(AutoCorrArrayTest, testSample)
TEST_ADD(AutoCorrArrayTest, testSerialize)
TEST_END(AutoCorrArrayTest)

#endif
