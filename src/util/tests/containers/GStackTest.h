#ifndef G_STACK_TEST_H
#define G_STACK_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <util/containers/DArray.h>
#include <util/containers/GStack.h>
#include <util/containers/PArrayIterator.h>

using namespace Util;

class GStackTest : public UnitTest
{

private:

   typedef int Data;

   const static int capacity = 10;
   DArray<Data>* arrayPtr;
   GStack<Data>* parrayPtr;
   int memory_;

   DArray<Data>& array()
   { return *arrayPtr; }

   GStack<Data>& stack()
   { return *parrayPtr; }

public:

   void setUp();
   void tearDown();

   void testPushPop();
   void testPushPopEmpty();

};

void GStackTest::setUp()
{
   memory_ = Memory::total();
   arrayPtr = new DArray<Data>;
   parrayPtr = new GStack<Data>;
   array().allocate(capacity);
   for (int i=0; i < capacity; i++) {
      array()[i] = (i+1)*10 + 1;
   }
}

void GStackTest::tearDown()
{}

void GStackTest::testPushPop()
{
   printMethod(TEST_FUNC);

   stack().reserve(2);
   stack().push(array()[8]); // 0
   stack().push(array()[4]); // 1
   TEST_ASSERT(stack().capacity() == 2);
   stack().push(array()[3]); // 2
   TEST_ASSERT(stack().capacity() == 4);
   stack().push(array()[5]); // 3
   TEST_ASSERT(stack().capacity() == 4);
   stack().push(array()[6]); // 4
   TEST_ASSERT(stack().capacity() == 8);
   stack().push(array()[9]); // 5

   TEST_ASSERT(stack().size() == 6);
   TEST_ASSERT(stack().pop() == array()[9]);
   TEST_ASSERT(stack().size() == 5);
   TEST_ASSERT(stack().pop() == array()[6]);
   TEST_ASSERT(stack().size() == 4);
   TEST_ASSERT(stack().peek() == array()[5]);
   TEST_ASSERT(stack().pop() == array()[5]);
   TEST_ASSERT(stack().size() == 3);
   TEST_ASSERT(stack().pop() == array()[3]);
   TEST_ASSERT(stack().size() == 2);
   TEST_ASSERT(stack().pop() == array()[4]);
   TEST_ASSERT(stack().size() == 1);
   TEST_ASSERT(stack().pop() == array()[8]);
   TEST_ASSERT(stack().size() == 0);
   TEST_ASSERT(stack().capacity() == 8);

   stack().clear();
   TEST_ASSERT(stack().size() == 0);
   TEST_ASSERT(stack().capacity() == 8);

   delete arrayPtr;
   delete parrayPtr;
   TEST_ASSERT(Memory::total() == memory_);
}

void GStackTest::testPushPopEmpty()
{
   printMethod(TEST_FUNC);
   {

      TEST_ASSERT(stack().size() == 0);
      TEST_ASSERT(stack().capacity() == 0);

      stack().push(array()[8]); // 0
      TEST_ASSERT(stack().capacity() == 64);
      TEST_ASSERT(stack().peek() == array()[8]);
      stack().push(array()[4]); // 1
      TEST_ASSERT(stack().capacity() == 64);
      TEST_ASSERT(stack().peek() == array()[4]);
      stack().push(array()[3]); // 2
      TEST_ASSERT(stack().peek() == array()[3]);
      stack().push(array()[5]); // 3
      stack().push(array()[6]); // 4
      stack().push(array()[9]); // 5
      TEST_ASSERT(stack().capacity() == 64);
      TEST_ASSERT(stack().isValid());

      TEST_ASSERT(stack().size() == 6);
      TEST_ASSERT(stack().pop() == array()[9]);
      TEST_ASSERT(stack().size() == 5);
      TEST_ASSERT(stack().pop() == array()[6]);
      TEST_ASSERT(stack().size() == 4);
      TEST_ASSERT(stack().peek() == array()[5]);
      TEST_ASSERT(stack().pop() == array()[5]);
      TEST_ASSERT(stack().size() == 3);
      TEST_ASSERT(stack().isValid());
      TEST_ASSERT(stack().pop() == array()[3]);
      TEST_ASSERT(stack().size() == 2);
      stack().push(array()[2]);
      TEST_ASSERT(stack().size() == 3);
      stack().push(array()[6]);
      TEST_ASSERT(stack().size() == 4);
      TEST_ASSERT(stack().pop() == array()[6]);
      TEST_ASSERT(stack().size() == 3);
      TEST_ASSERT(stack().pop() == array()[2]);
      TEST_ASSERT(stack().size() == 2);
      TEST_ASSERT(stack().pop() == array()[4]);
      TEST_ASSERT(stack().size() == 1);
      TEST_ASSERT(stack().pop() == array()[8]);
      TEST_ASSERT(stack().size() == 0);

      delete arrayPtr;
      delete parrayPtr;
   }
   TEST_ASSERT(Memory::total() == memory_);
}

TEST_BEGIN(GStackTest)
TEST_ADD(GStackTest, testPushPop)
TEST_ADD(GStackTest, testPushPopEmpty)
TEST_END(GStackTest)

#endif
