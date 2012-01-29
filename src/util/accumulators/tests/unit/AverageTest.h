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

class AverageTest : public ParamFileTest<Average>
{

public:

   AverageTest()
   { } 

   void setUp() {
      openFile("in/Average"); 
      object().readParam(file());
      closeFile();
   }
      
   void readData() 
   {
      int i, n;
      double x;
      std::ifstream datafile("in/data");
      datafile >> n;
      for (i = 0; i < n; ++i) {
         datafile >> x;
         object().sample(x);
      }
      datafile.close();
   }

   void testReadParam() 
   {
      printMethod(TEST_FUNC);

      openFile("in/Average"); 
      object().readParam(file());
      closeFile();

      printEndl();      
      object().writeParam(std::cout);
   }

   void testSample() 
   {
      printMethod(TEST_FUNC);

      readData();

      printEndl();      
      object().output(std::cout);
   }

   void testSerialize() 
   {
      printMethod(TEST_FUNC);
      printEndl();

      readData();

      MemoryOArchive u;
      int size = memorySize(object());
      u.allocate(size);

      std::cout << size << std::endl;

      u << object();
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
