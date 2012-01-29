#ifndef IO_UTIL_TEST_H
#define IO_UTIL_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <util/util/ioUtil.h>
#include <util/global.h>

using namespace Util;

class ioUtilTest : public UnitTest 
{

public:

   void setUp()
   {};

   void tearDown()
   {};

   void testGetNextLine() 
   {
      printMethod(TEST_FUNC);

      std::string line;
      std::istream in("in/GetNextLine");
      getNextLine(in, line);
      std::cout << line << std::endl;
   }

};

TEST_BEGIN(ioUtilTest)
TEST_ADD(ioUtilTest, testGetNextLine)
TEST_END(ioUtilTest)

#endif
