#ifndef IO_UTIL_TEST_H
#define IO_UTIL_TEST_H

#include <util/util/ioUtil.h>
#include <util/global.h>

#ifdef UTIL_MPI
#ifndef TEST_MPI
#define TEST_MPI
#endif
#endif

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

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
      std::ifstream in("in/GetNextLine");
      getNextLine(in, line);
      std::cout << line << std::endl;
   }

   void testRStrip() 
   {
      printMethod(TEST_FUNC);

      std::string full("  This has white space   ");
      std::string lean("  This has white space");
      int len = rStrip(full);
      TEST_ASSERT(len == 22);
      TEST_ASSERT(full == lean);
   }

};

TEST_BEGIN(ioUtilTest)
//TEST_ADD(ioUtilTest, testGetNextLine)
TEST_ADD(ioUtilTest, testRStrip)
TEST_END(ioUtilTest)

#endif
