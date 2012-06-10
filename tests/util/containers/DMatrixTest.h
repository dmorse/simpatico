#ifndef DMATRIX_TEST_H
#define DMATRIX_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <util/containers/DMatrix.h>
#include <util/containers/Matrix.h>

using namespace Util;

class DMatrixTest : public UnitTest 
{

public:

   void setUp() {}

   void tearDown() {}
  
   void testConstructor();

   void testAllocate();

   void testSubscript();

   void testCopyConstructor();

   void testAssignment();

   void testBaseClassReference();

};


void DMatrixTest::testConstructor()
{
   printMethod(TEST_FUNC);
   DMatrix<int> v;
   TEST_ASSERT(v.capacity1() == 0 );
   TEST_ASSERT(v.capacity2() == 0 );
   TEST_ASSERT(!v.isAllocated() );
} 

void DMatrixTest::testAllocate()
{
   printMethod(TEST_FUNC);
   DMatrix<int> v;
   v.allocate(2,3);
   TEST_ASSERT(v.capacity1() == 2 );
   TEST_ASSERT(v.capacity2() == 3 );
   TEST_ASSERT(v.isAllocated() );
} 

void DMatrixTest::testSubscript()
{
   printMethod(TEST_FUNC);
   DMatrix<int> v;
   v.allocate(2,2);
   v(0,0) = 3;
   v(1,0) = 4;
   v(0,1) = 5;
   v(1,1) = 6;
   TEST_ASSERT(v(0,0) == 3 );
   TEST_ASSERT(v(1,0) == 4 );
   TEST_ASSERT(v(0,1) == 5 );
   TEST_ASSERT(v(1,1) == 6 );
} 

void DMatrixTest::testCopyConstructor()
{
   printMethod(TEST_FUNC);
   DMatrix<int> v;
   TEST_ASSERT(v.capacity1() == 0 );
   TEST_ASSERT(v.capacity2() == 0 );
   TEST_ASSERT(!v.isAllocated() );

   v.allocate(2,2);
   TEST_ASSERT(v.capacity1() == 2 );
   TEST_ASSERT(v.capacity2() == 2 );
   TEST_ASSERT(v.isAllocated() );
   v(0,0) = 3;
   v(1,0) = 4;
   v(0,1) = 5;
   v(1,1) = 6;

   DMatrix<int> u(v);
   TEST_ASSERT(u.capacity1() == 2 );
   TEST_ASSERT(u.capacity2() == 2 );
   TEST_ASSERT(u.isAllocated() );
   TEST_ASSERT(v(0,0) == 3 );
   TEST_ASSERT(v(1,0) == 4 );
   TEST_ASSERT(v(0,1) == 5 );
   TEST_ASSERT(v(1,1) == 6 );
   TEST_ASSERT(u(0,0) == 3 );
   TEST_ASSERT(u(1,0) == 4 );
   TEST_ASSERT(u(0,1) == 5 );
   TEST_ASSERT(u(1,1) == 6 );
} 

void DMatrixTest::testAssignment()
{
   printMethod(TEST_FUNC);

   DMatrix<int> v;
   v.allocate(2,2);
   TEST_ASSERT(v.capacity1() == 2);
   TEST_ASSERT(v.capacity2() == 2);
   TEST_ASSERT(v.isAllocated() );

   DMatrix<int> u;
   u.allocate(2,2);
   TEST_ASSERT(v.capacity1() == 2);
   TEST_ASSERT(v.capacity2() == 2);
   TEST_ASSERT(v.isAllocated() );

   v(0,0) = 3;
   v(1,0) = 4;
   v(0,1) = 5;
   v(1,1) = 6;

   u  = v;

   TEST_ASSERT(u.capacity1() == 2 );
   TEST_ASSERT(u.capacity2() == 2 );
   TEST_ASSERT(u.isAllocated() );

   TEST_ASSERT(v(0,0) == 3 );
   TEST_ASSERT(v(1,0) == 4 );
   TEST_ASSERT(v(0,1) == 5 );
   TEST_ASSERT(v(1,1) == 6 );

   TEST_ASSERT(u(0,0) == 3 );
   TEST_ASSERT(u(1,0) == 4 );
   TEST_ASSERT(u(0,1) == 5 );
   TEST_ASSERT(u(1,1) == 6 );
} 

void DMatrixTest::testBaseClassReference()
{
   printMethod(TEST_FUNC);
   DMatrix<int> v;
   v.allocate(2,2);
   v(0,0) = 3;
   v(1,0) = 4;
   v(0,1) = 5;
   v(1,1) = 6;

   Matrix<int>& u = v;

   TEST_ASSERT(u.capacity1() == 2 );
   TEST_ASSERT(u.capacity2() == 2 );
   TEST_ASSERT(u(0,0) == 3 );
   TEST_ASSERT(u(1,0) == 4 );
   TEST_ASSERT(u(0,1) == 5 );
   TEST_ASSERT(u(1,1) == 6 );

} 

TEST_BEGIN(DMatrixTest)
TEST_ADD(DMatrixTest, testConstructor)
TEST_ADD(DMatrixTest, testAllocate)
TEST_ADD(DMatrixTest, testSubscript)
TEST_ADD(DMatrixTest, testCopyConstructor)
TEST_ADD(DMatrixTest, testAssignment)
TEST_ADD(DMatrixTest, testBaseClassReference)
TEST_END(DMatrixTest)

#endif
