#ifndef LABEL_TEST_H
#define LABEL_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <util/param/Label.h>

#include <iostream>
#include <fstream>

using namespace Util;

class LabelTest : public UnitTest 
{

public:

   void setUp()
   { }

   void tearDown()
   { }

   void testLabelConstructor() 
   {
      printMethod(TEST_FUNC);
      Label label("MyLabel");
   }


   void testExtracter() 
   {
      printMethod(TEST_FUNC);
      Label label("MyLabel");
      std::ifstream in;
      openInputFile("in/Label", in);
      in >> label;
      in >> label;
      in.close();
   }


   void testInserter() 
   {
      printMethod(TEST_FUNC);
      printEndl();
      Label label("MyLabel");
      std::cout << label;
      std::cout << "\n"; 
   }

};


TEST_BEGIN(LabelTest)
TEST_ADD(LabelTest, testLabelConstructor)
TEST_ADD(LabelTest, testExtracter)
TEST_ADD(LabelTest, testInserter)
TEST_END(LabelTest)

#endif
