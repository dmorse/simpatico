#ifndef AUTOCORR_TEST_H
#define AUTOCORR_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <util/accumulators/AutoCorr.h>
#include <util/archives/MemoryOArchive.h>
#include <util/archives/MemoryIArchive.h>
#include <util/archives/MemoryCounter.h>

#include <iostream>
#include <fstream>

using namespace Util;

class AutoCorrTest : public UnitTest
{

   AutoCorr<double, double> accumulator_;

public:

   AutoCorrTest() 
   {}

   void setUp() 
   {
      std::ifstream paramFile; 
      openInputFile("in/AutoCorr", paramFile); 
      accumulator_.readParam(paramFile);
      paramFile.close();
   }

   void tearDown() 
   {}

   void readData() 
   {
      int i, n;
      double x;
      std::ifstream dataFile; 
      openInputFile("in/data", dataFile); 
      dataFile >> n;
      for (i = 0; i < n; ++i) {
         dataFile >> x;
         accumulator_.sample(x);
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

   #if 0
   void testPack() 
   {
      printMethod(TEST_FUNC);
      printEndl();

      readData();

      int size = accumulator_.packedSize();
      std::cout << size << std::endl;

      char* begin  = new char[size];
      char* end    = begin + size;
      char* cursor = begin;
      accumulator_.pack(cursor, end);
      TEST_ASSERT(cursor == end);

      AutoCorr<double, double> clone;
      clone.setParam(32);
      cursor = begin;
      clone.unpack(cursor, end);

      clone.output(std::cout);

   }
   #endif

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

      AutoCorr<double, double> clone;
      v & clone;

      clone.output(std::cout);
   }

};

TEST_BEGIN(AutoCorrTest)
TEST_ADD(AutoCorrTest, testReadParam)
TEST_ADD(AutoCorrTest, testSample)
//TEST_ADD(AutoCorrTest, testPack)
TEST_ADD(AutoCorrTest, testSerialize)
TEST_END(AutoCorrTest)

#endif
