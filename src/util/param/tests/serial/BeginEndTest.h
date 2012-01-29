#ifndef BEGIN_END_TEST_H
#define BEGIN_END_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <util/param/Begin.h>
#include <util/param/End.h>

#include <iostream>
#include <fstream>

using namespace Util;

class BeginEndTest : public UnitTest 
{

public:

   void setUp()
   {}

   void tearDown()
   {}

   void testBeginConstructor() {
      printMethod(TEST_FUNC);
      Begin* param;
      param = new Begin("ClassName");
      delete param;
   }

   void testBeginWrite() {
      printMethod(TEST_FUNC);
      Begin* param;
      param = new Begin("ClassName");
      if (verbose() > 0) {
         param->writeParam(std::cout);
      }
      delete param;
   }

   void testBeginRead() {
      printMethod(TEST_FUNC);
      Begin *param;
      param = new Begin("ClassName");
      std::ifstream in;
      openInputFile("in/Begin", in);
      param->readParam(in);
      if (verbose() > 0) {
         param->writeParam(std::cout);
      }
      delete param;
   }

   void testEndWrite() {
      printMethod(TEST_FUNC);
      End* param;
      param = new End();
      if (verbose() > 0) {
         param->writeParam(std::cout);
      }
      delete param;
   }

   void testEndRead() {
      printMethod(TEST_FUNC);
      End *param;
      param = new End();
      std::ifstream in;
      openInputFile("in/End", in);
      param->readParam(in);
      if (verbose() > 0) {
         param->writeParam(std::cout);
      }
      delete param;
   }

};

TEST_BEGIN(BeginEndTest)
TEST_ADD(BeginEndTest, testBeginConstructor)
TEST_ADD(BeginEndTest, testBeginWrite)
TEST_ADD(BeginEndTest, testBeginRead)
TEST_ADD(BeginEndTest, testEndWrite)
TEST_ADD(BeginEndTest, testEndRead)
TEST_END(BeginEndTest)

#endif
