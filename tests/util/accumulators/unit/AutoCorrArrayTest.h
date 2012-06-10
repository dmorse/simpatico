#ifndef AUTOCORR_ARRAY_TEST_H
#define AUTOCORR_ARRAY_TEST_H

#include <test/UnitTest.h>
#include <test/ParamFileTest.h>
#include <test/UnitTestRunner.h>

#include <util/accumulators/AutoCorrArray.h>
#include <util/archives/MemoryOArchive.h>
#include <util/archives/MemoryIArchive.h>
#include <util/archives/MemoryCounter.h>

#include <iostream>
#include <fstream>

using namespace Util;

class AutoCorrArrayTest : public ParamFileTest< AutoCorrArray<double, double> >
{

public:

   //typedef AutoCorrArray<double, double> Object;

   AutoCorrArrayTest() 
   {}

   void setUp() 
   { 
      openFile("in/AutoCorrArray"); 
      object().readParam(file());
      closeFile();
      object().setNEnsemble(4);
   }

   void tearDown() 
   {}

   void readData() 
   {
      DArray<double> data;
      int i, j, m, n, p;

      std::ifstream datafile("in/data");
      datafile >> m;
      n = object().nEnsemble();
      data.allocate(n);
      p = m / n;
      for (i = 0; i < p; ++i) {
         for (j = 0; j < n; ++j) {
           datafile >> data[j];
         }
         object().sample(data);
      }
      datafile.close();
   }

   void testReadParam() 
   {
      printMethod(TEST_FUNC);

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

      int size = memorySize(object());

      MemoryOArchive u;
      u.allocate(size);

      std::cout << size << std::endl;

      u << object();
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
