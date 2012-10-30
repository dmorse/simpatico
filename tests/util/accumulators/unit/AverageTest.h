#ifndef AVERAGE_TEST_H
#define AVERAGE_TEST_H

#include <test/UnitTest.h>
#include <test/ParamFileTest.h>
#include <test/UnitTestRunner.h>

#include <util/accumulators/Average.h>
#include <util/archives/MemoryOArchive.h>
#include <util/archives/MemoryIArchive.h>
#include <util/archives/MemoryCounter.h>

using namespace Util;

class AverageTest : public UnitTest
{

   Average accumulator_;

public:

   AverageTest()
   { } 

   void setUp() {
      std::ifstream paramFile;
      openInputFile("in/Average", paramFile); 
      accumulator_.readParam(paramFile);
      paramFile.close();
   }
      
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

   void testSerialize() 
   {
      printMethod(TEST_FUNC);
      printEndl();

      readData();

      MemoryOArchive u;
      int size = memorySize(accumulator_);
      u.allocate(size);

      std::cout << size << std::endl;

      u << accumulator_;
      TEST_ASSERT(u.cursor() == u.begin() + size);

      MemoryIArchive v;
      v = u;

      Average clone;
      v >> clone;
      TEST_ASSERT(v.cursor() == u.begin() + size);

      clone.output(std::cout);
   }


};

TEST_BEGIN(AverageTest)
TEST_ADD(AverageTest, testReadParam)
TEST_ADD(AverageTest, testSample)
TEST_ADD(AverageTest, testSerialize)
TEST_END(AverageTest)

#endif
