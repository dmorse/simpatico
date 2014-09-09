#ifndef XML_ATTRIBUTE_TEST_H
#define XML_ATTRIBUTE_TEST_H

#include <util/misc/Xml.h>
#include <util/global.h>

#ifdef UTIL_MPI
#ifndef TEST_MPI
#define TEST_MPI
#endif
#endif

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

using namespace Util;

class XmlTest : public UnitTest 
{

public:

   void setUp()
   {};

   void tearDown()
   {};

   void testMatch1() 
   {
      printMethod(TEST_FUNC);

      XmlAttribute parser;
      std::string string("  this =\"35.5\"");
      bool result = parser.match(string, 0);
      TEST_ASSERT(result == true);
      //std::cout << "\n" << parser.label();
      //std::cout << "\n" << parser.value().str();
      double value;
      parser.value() >>  value;
      TEST_ASSERT(parser.label() == std::string("this"));
      TEST_ASSERT(value == 35.5);
   }

   void testMatch2() 
   {
      printMethod(TEST_FUNC);

      ParserString parent;
      XmlAttribute parser;
      std::string string("  this =\"35.5\"u");
      parent.setString(string, 0);
      bool result = parser.match(parent);
      TEST_ASSERT(result == true);
      //std::cout << "\n" << parser.label();
      //std::cout << "\n" << parser.value().str();
      double value;
      parser.value() >>  value;
      TEST_ASSERT(parser.label() == std::string("this"));
      TEST_ASSERT(value == 35.5);
      TEST_ASSERT(parser.cursor() == parent.cursor());
      TEST_ASSERT(parent.c() == 'u');
      TEST_ASSERT(parent.cursor() == 14);
   }

   void testXmlStartTag() 
   {
      printMethod(TEST_FUNC);

      XmlStartTag tag;
      std::string string("<Label this =\"35.5\" >");
      bool result = tag.matchLabel(string, 0);
      TEST_ASSERT(result);
   }

};

TEST_BEGIN(XmlTest)
TEST_ADD(XmlTest, testMatch1)
TEST_ADD(XmlTest, testMatch2)
TEST_ADD(XmlTest, testXmlStartTag)
TEST_END(XmlTest)

#endif
