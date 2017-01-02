#ifndef FLAG_SET_TEST_H
#define FLAG_SET_TEST_H

#include <util/misc/FlagSet.h>
#include <util/global.h>

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <string>

using namespace Util;

class FlagSetTest : public UnitTest 
{

public:

   void setUp()
   {};

   void tearDown()
   {};

   void testConstructor() 
   {
      printMethod(TEST_FUNC);

      std::string allowed = "abcdef";
      FlagSet flags(allowed);
      TEST_ASSERT(flags.allowed() == allowed);
      TEST_ASSERT(flags.isActive('a') == false);
      TEST_ASSERT(flags.isActive('b') == false);
      TEST_ASSERT(flags.isActive('c') == false);
      TEST_ASSERT(flags.isActive('d') == false);
      TEST_ASSERT(flags.isActive('e') == false);
      TEST_ASSERT(flags.isActive('f') == false);
   }

   void testSetAllowed() 
   {
      printMethod(TEST_FUNC);

      FlagSet flags;
      std::string allowed = "abcdef";
      flags.setAllowed(allowed);
      TEST_ASSERT(flags.allowed() == allowed);
      TEST_ASSERT(flags.isActive('a') == false);
      TEST_ASSERT(flags.isActive('b') == false);
      TEST_ASSERT(flags.isActive('c') == false);
      TEST_ASSERT(flags.isActive('d') == false);
      TEST_ASSERT(flags.isActive('e') == false);
      TEST_ASSERT(flags.isActive('f') == false);
   }

   void testReadOrdered1() 
   {
      printMethod(TEST_FUNC);

      FlagSet flags;
      std::string allowed = "abcdef";
      flags.setAllowed(allowed);
      TEST_ASSERT(flags.allowed() == allowed);
      std::string actual = "bde";
      flags.setActualOrdered(actual);
      TEST_ASSERT(flags.actual() == actual);
      TEST_ASSERT(flags.isActive('a') == false);
      TEST_ASSERT(flags.isActive('b') == true);
      TEST_ASSERT(flags.isActive('c') == false);
      TEST_ASSERT(flags.isActive('d') == true);
      TEST_ASSERT(flags.isActive('e') == true);
      TEST_ASSERT(flags.isActive('f') == false);
   }

   void testReadOrdered2() 
   {
      printMethod(TEST_FUNC);

      FlagSet flags;
      std::string allowed = "abcdef";
      flags.setAllowed(allowed);
      TEST_ASSERT(flags.allowed() == allowed);
      std::string actual = "acf";
      flags.setActualOrdered(actual);
      TEST_ASSERT(flags.actual() == actual);
      TEST_ASSERT(flags.isActive('a') == true);
      TEST_ASSERT(flags.isActive('b') == false);
      TEST_ASSERT(flags.isActive('c') == true);
      TEST_ASSERT(flags.isActive('d') == false);
      TEST_ASSERT(flags.isActive('e') == false);
      TEST_ASSERT(flags.isActive('f') == true);
   }

   void testReadOrdered3() 
   {
      printMethod(TEST_FUNC);

      FlagSet flags;
      std::string allowed = "abcdef";
      flags.setAllowed(allowed);
      TEST_ASSERT(flags.allowed() == allowed);
      bool success = false;
      std::string actual = "afc";
      try {
         flags.setActualOrdered(actual);
      } catch (Exception e) {
         success = true;
      }
      TEST_ASSERT(success);
   }

   void testReadOrdered4() 
   {
      printMethod(TEST_FUNC);

      FlagSet flags;
      std::string allowed = "abcdef";
      flags.setAllowed(allowed);
      TEST_ASSERT(flags.allowed() == allowed);
      bool success = false;
      std::string actual = "ag";
      try {
         flags.setActualOrdered(actual);
      } catch (Exception e) {
         success = true;
      }
      TEST_ASSERT(success);
   }

};

TEST_BEGIN(FlagSetTest)
TEST_ADD(FlagSetTest, testConstructor)
TEST_ADD(FlagSetTest, testSetAllowed)
TEST_ADD(FlagSetTest, testReadOrdered1)
TEST_ADD(FlagSetTest, testReadOrdered2)
TEST_ADD(FlagSetTest, testReadOrdered3)
TEST_ADD(FlagSetTest, testReadOrdered4)
TEST_END(FlagSetTest)

#endif
