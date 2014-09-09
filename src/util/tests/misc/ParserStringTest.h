#ifndef PARSER_STRING_TEST_H
#define PARSER_STRING_TEST_H

#include <util/misc/ParserString.h>
#include <util/global.h>

#ifdef UTIL_MPI
#ifndef TEST_MPI
#define TEST_MPI
#endif
#endif

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

using namespace Util;

class ParserStringTest : public UnitTest 
{

public:

   void setUp()
   {};

   void tearDown()
   {};

   void testSetString() 
   {
      printMethod(TEST_FUNC);

      ParserString parser;
      std::string string("this is the string");
      parser.setString("this is the string", 0);
      TEST_ASSERT(parser.string() == string);
      TEST_ASSERT(parser.cursor() == 0);
      TEST_ASSERT(parser.c() == 't');
   }

   void testMove() 
   {
      printMethod(TEST_FUNC);

      ParserString parser;
      std::string string("   this is the string");
      parser.setString(string, 0);
      TEST_ASSERT(parser.string() == string);
      TEST_ASSERT(parser.cursor() == 0);
      TEST_ASSERT(parser.c() == ' ');
      parser.skip();
      TEST_ASSERT(parser.cursor() == 3);
      TEST_ASSERT(parser.c() == 't');
      parser.next();
      TEST_ASSERT(parser.cursor() == 4);
      TEST_ASSERT(parser.c() == 'h');
      while (!parser.isEnd()) {
         parser.next();
      }
      TEST_ASSERT(parser.isEnd());
      TEST_ASSERT(parser.c() == '\0');
      TEST_ASSERT(parser.cursor() == 21);
      
   }

};

TEST_BEGIN(ParserStringTest)
TEST_ADD(ParserStringTest, testSetString)
TEST_ADD(ParserStringTest, testMove)
TEST_END(ParserStringTest)

#endif
