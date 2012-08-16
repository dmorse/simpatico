#ifndef NULLABLE_TEST_H
#define NULLABLE_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>
#include <util/util/Nullable.h>
#include <util/global.h>

using namespace Util;

class NullableTest : public UnitTest 
{

public:

   void setUp()
   {};

   void tearDown()
   {};

   void testDefaultConstructor() 
   {
      printMethod(TEST_FUNC);

      Nullable<int> nullable;
      TEST_ASSERT(nullable.isNull());
   }

   void testConstructorFromValue() 
   {
      printMethod(TEST_FUNC);

      Nullable<int> nullable(3);
      TEST_ASSERT(!nullable.isNull());
      TEST_ASSERT(nullable.value() == 3);
      nullable.setNull();
      TEST_ASSERT(nullable.isNull());
   }

   void testAssignment() 
   {
      printMethod(TEST_FUNC);

      Nullable<int> nullable1(3);
      TEST_ASSERT(!nullable1.isNull());
      TEST_ASSERT(nullable1.value() == 3);

      Nullable<int> nullable2;
      TEST_ASSERT(nullable2.isNull());
        
      nullable2 = nullable1;
      TEST_ASSERT(!nullable2.isNull());
      TEST_ASSERT(nullable2.value() == 3);
   }

   void testAssignmentFromValue() 
   {
      printMethod(TEST_FUNC);

      Nullable<int> nullable;
      TEST_ASSERT(nullable.isNull());

      nullable = 3;
      TEST_ASSERT(!nullable.isNull());
      TEST_ASSERT(nullable.value() == 3);

      nullable.setNull();
      TEST_ASSERT(nullable.isNull());
   }

};

TEST_BEGIN(NullableTest)
TEST_ADD(NullableTest, testDefaultConstructor)
TEST_ADD(NullableTest, testConstructorFromValue)
TEST_ADD(NullableTest, testAssignment)
TEST_ADD(NullableTest, testAssignmentFromValue)
TEST_END(NullableTest)

#endif
