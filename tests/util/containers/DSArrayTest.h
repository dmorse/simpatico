#ifndef DS_ARRAY_TEST_H
#define DS_ARRAY_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <util/containers/DSArray.h>
#include <util/containers/Array.h>

using namespace Util;

class DSArrayTest : public UnitTest 
{

public:

   void setUp(){};
   void tearDown(){};
   void testAllocate();
   void testConstructor();
   void testSubscript();
   void testCopyConstructor();

};


void DSArrayTest::testConstructor()
{
   printMethod(TEST_FUNC);
   DSArray<int> v;
   TEST_ASSERT(v.capacity() == 0);
   TEST_ASSERT(!v.isAllocated() );
} 

void DSArrayTest::testAllocate()
{
   printMethod(TEST_FUNC);
   DSArray<int> v;
   v.allocate(3);
   TEST_ASSERT(v.capacity() == 3 );
   TEST_ASSERT(v.isAllocated() );
   TEST_ASSERT(v.size() == 0 );
}

void DSArrayTest::testSubscript()
{
   printMethod(TEST_FUNC);
   DSArray<int> v;
   v.allocate(3);
   v.append(3);
   v.append(4);
   v.append(5);
   TEST_ASSERT(v.size() == 3);
   TEST_ASSERT(v[0] == 3);
   TEST_ASSERT(v[1] == 4);
   TEST_ASSERT(v[2] == 5);
} 

void DSArrayTest::testCopyConstructor()
{
   printMethod(TEST_FUNC);
   DSArray<int> v;
   TEST_ASSERT(v.capacity() == 0 );
   TEST_ASSERT(!v.isAllocated() );

   v.allocate(3);
   TEST_ASSERT(v.capacity() == 3 );
   TEST_ASSERT(v.isAllocated() );
   v.append(3);
   v.append(4);
   v.append(5);

   DSArray<int> u(v);
   TEST_ASSERT(u.capacity() == 3 );
   TEST_ASSERT(u.isAllocated() );
   TEST_ASSERT(u.size() == 3);
   TEST_ASSERT(v.size() == 3);
   TEST_ASSERT(v[0] == 3 );
   TEST_ASSERT(v[1] == 4 );
   TEST_ASSERT(v[2] == 5 );
   TEST_ASSERT(u[0] == 3 );
   TEST_ASSERT(u[1] == 4 );
   TEST_ASSERT(u[2] == 5 );
}


TEST_BEGIN(DSArrayTest)
TEST_ADD(DSArrayTest, testAllocate)
TEST_ADD(DSArrayTest, testConstructor)
TEST_ADD(DSArrayTest, testSubscript)
TEST_ADD(DSArrayTest, testCopyConstructor)
TEST_END(DSArrayTest)

#endif
