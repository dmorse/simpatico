#ifndef AUTOCORR_TEST_H
#define AUTOCORR_TEST_H

#include <test/UnitTest.h>
#include <test/ParamFileTest.h>
#include <test/UnitTestRunner.h>

#include <util/accumulators/AutoCorr.h>
#include <util/archives/MemoryOArchive.h>
#include <util/archives/MemoryIArchive.h>
#include <util/archives/MemoryCounter.h>

#include <iostream>
#include <fstream>

using namespace Util;

class AutoCorrTest : public ParamFileTest< AutoCorr<double, double> >
{

public:

   //typedef AutoCorr<double, double> Object;

   AutoCorrTest() 
   {}

   void setUp() 
   { 
      openFile("in/AutoCorr"); 
      object().readParam(file());
      closeFile();
   }

   void tearDown() 
   {}

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

   void testPack() 
   {
      printMethod(TEST_FUNC);
      printEndl();

      readData();

      int size = object().packedSize();
      std::cout << size << std::endl;

      char* begin  = new char[size];
      char* end    = begin + size;
      char* cursor = begin;
      object().pack(cursor, end);
      TEST_ASSERT(cursor == end);

      AutoCorr<double, double> clone;
      clone.setParam(32);
      cursor = begin;
      clone.unpack(cursor, end);

      clone.output(std::cout);

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

      AutoCorr<double, double> clone;
      v & clone;

      clone.output(std::cout);
   }

};

TEST_BEGIN(AutoCorrTest)
TEST_ADD(AutoCorrTest, testReadParam)
TEST_ADD(AutoCorrTest, testSample)
TEST_ADD(AutoCorrTest, testPack)
TEST_ADD(AutoCorrTest, testSerialize)
TEST_END(AutoCorrTest)

#endif
