#ifndef GRID_ARRAY_TEST_H
#define GRID_ARRAY_TEST_H

#include <util/containers/GridArray.h>
#include <util/space/IntVector.h>

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

using namespace Util;

class GridArrayTest : public UnitTest 
{

public:

   void setUp() {}

   void tearDown() {}
  
   void testConstructor();

   void testAllocate();

   void testSubscript();

   #if 0
   void testCopyConstructor();

   void testAssignment();

   void testSerializeFile1();

   void testSerializeFile2();
   #endif

};


void GridArrayTest::testConstructor()
{
   printMethod(TEST_FUNC);
   GridArray<int> v;
   TEST_ASSERT(!v.isAllocated() );
} 

void GridArrayTest::testAllocate()
{
   printMethod(TEST_FUNC);
   GridArray<int> v;
   TEST_ASSERT(!v.isAllocated());

   IntVector dimensions;
   dimensions[0] = 4;
   dimensions[1] = 3;
   dimensions[2] = 2;
   v.allocate(dimensions);

   TEST_ASSERT(v.size() == 24);
   TEST_ASSERT(v.dimension(0) == 4);
   TEST_ASSERT(v.dimension(1) == 3);
   TEST_ASSERT(v.dimension(2) == 2);
   TEST_ASSERT(v.dimensions() == dimensions);
   TEST_ASSERT(v.isAllocated());

   IntVector position;
   position[0] = 2;
   position[1] = 1;
   position[2] = 1;
   TEST_ASSERT(v.rank(position) == 15);
   TEST_ASSERT(v.position(15) == position);

   position[0] = 3;
   position[1] = 2;
   position[2] = 1;
   TEST_ASSERT(v.rank(position) == 23);
   TEST_ASSERT(v.position(23) == position);

} 

void GridArrayTest::testSubscript()
{
   printMethod(TEST_FUNC);
   GridArray<int> v;

   IntVector dimensions;
   dimensions[0] = 4;
   dimensions[1] = 3;
   dimensions[2] = 2;
   v.allocate(dimensions);

   TEST_ASSERT(v.size() == 24);
   TEST_ASSERT(v.dimension(0) == 4);
   TEST_ASSERT(v.dimension(1) == 3);
   TEST_ASSERT(v.dimension(2) == 2);
   TEST_ASSERT(v.dimensions() == dimensions);
   TEST_ASSERT(v.isAllocated());

   IntVector position;
   position[0] = 2;
   position[1] = 1;
   position[2] = 1;
   v(position) = 38;
   TEST_ASSERT(v.rank(position) == 15);
   TEST_ASSERT(v.position(15) == position);
   TEST_ASSERT(v[15] == 38);
   TEST_ASSERT(v(position) == 38);
} 

#if 0
void GridArrayTest::testCopyConstructor()
{
   printMethod(TEST_FUNC);
   GridArray<int> v;
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

   GridArray<int> u(v);
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

void GridArrayTest::testAssignment()
{
   printMethod(TEST_FUNC);

   GridArray<int> v;
   v.allocate(2,2);
   TEST_ASSERT(v.capacity1() == 2);
   TEST_ASSERT(v.capacity2() == 2);
   TEST_ASSERT(v.isAllocated() );

   GridArray<int> u;
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

void GridArrayTest::testSerializeFile1()
{
   printMethod(TEST_FUNC);
   GridArray<int> v;
   v.allocate(2,2);
   v(0,0) = 3;
   v(1,0) = 4;
   v(0,1) = 5;
   v(1,1) = 6;
   int i1 = 13;
   int i2;

   TEST_ASSERT(v.capacity1() == 2 );
   TEST_ASSERT(v.capacity2() == 2 );
   TEST_ASSERT(v(0,0) == 3 );
   TEST_ASSERT(v(1,0) == 4 );
   TEST_ASSERT(v(0,1) == 5 );
   TEST_ASSERT(v(1,1) == 6 );

   BinaryFileOArchive oArchive;
   openOutputFile("binary", oArchive.file());
   oArchive << v;
   oArchive << i1;
   oArchive.file().close();

   // Show that v is unchanged by packing
   TEST_ASSERT(v.capacity1() == 2 );
   TEST_ASSERT(v.capacity2() == 2 );
   TEST_ASSERT(v(0,0) == 3 );
   TEST_ASSERT(v(1,0) == 4 );
   TEST_ASSERT(v(0,1) == 5 );
   TEST_ASSERT(v(1,1) == 6 );

   GridArray<int> u;
   u.allocate(2, 2);

   BinaryFileIArchive iArchive;
   openInputFile("binary", iArchive.file());
   iArchive >> u;
   iArchive >> i2;
   iArchive.file().close();

   TEST_ASSERT(u.capacity1() == 2);
   TEST_ASSERT(u.capacity2() == 2);
   TEST_ASSERT(u(0,0) == 3);
   TEST_ASSERT(u(1,0) == 4);
   TEST_ASSERT(u(0,1) == 5);
   TEST_ASSERT(u(1,1) == 6);
   TEST_ASSERT(i2 == i1);
   TEST_ASSERT(i2 == 13);

   #if 0
   // Clear values of u and i2
   for (int i=0; i < capacity; i++ ) {
      real(u[i]) = 0.0;
      imag(u[i]) = 0.0;
   }
   i2 = 0;

   // Reload into u and i2
   openInputFile("binary", iArchive.file());
   iArchive >> u;
   iArchive >> i2;

   TEST_ASSERT(imag(u[0]) == 10.1);
   TEST_ASSERT(real(u[1]) == 20.0);
   TEST_ASSERT(imag(u[2]) == 30.1);
   TEST_ASSERT(i2 == 13);
   TEST_ASSERT(u.capacity() == 3);
   #endif

} 

void GridArrayTest::testSerializeFile2()
{
   printMethod(TEST_FUNC);
   GridArray<int> v;
   v.allocate(2,2);
   v(0,0) = 3;
   v(1,0) = 4;
   v(0,1) = 5;
   v(1,1) = 6;
   int i1 = 13;
   int i2;

   TEST_ASSERT(v.capacity1() == 2 );
   TEST_ASSERT(v.capacity2() == 2 );
   TEST_ASSERT(v(0,0) == 3 );
   TEST_ASSERT(v(1,0) == 4 );
   TEST_ASSERT(v(0,1) == 5 );
   TEST_ASSERT(v(1,1) == 6 );

   BinaryFileOArchive oArchive;
   openOutputFile("binary", oArchive.file());
   oArchive << v;
   oArchive << i1;
   oArchive.file().close();

   // Show that v is unchanged by packing
   TEST_ASSERT(v.capacity1() == 2 );
   TEST_ASSERT(v.capacity2() == 2 );
   TEST_ASSERT(v(0,0) == 3 );
   TEST_ASSERT(v(1,0) == 4 );
   TEST_ASSERT(v(0,1) == 5 );
   TEST_ASSERT(v(1,1) == 6 );

   GridArray<int> u;

   //u.allocate(2, 2);  -> Note allocation, different from previous

   BinaryFileIArchive iArchive;
   openInputFile("binary", iArchive.file());
   iArchive >> u;
   iArchive >> i2;
   iArchive.file().close();

   TEST_ASSERT(u.capacity1() == 2);
   TEST_ASSERT(u.capacity2() == 2);
   TEST_ASSERT(u(0,0) == 3);
   TEST_ASSERT(u(1,0) == 4);
   TEST_ASSERT(u(0,1) == 5);
   TEST_ASSERT(u(1,1) == 6);
   TEST_ASSERT(i2 == i1);
   TEST_ASSERT(i2 == 13);

   #if 0
   // Clear values of u and i2
   for (int i=0; i < capacity; i++ ) {
      real(u[i]) = 0.0;
      imag(u[i]) = 0.0;
   }
   i2 = 0;

   // Reload into u and i2
   openInputFile("binary", iArchive.file());
   iArchive >> u;
   iArchive >> i2;

   TEST_ASSERT(imag(u[0]) == 10.1);
   TEST_ASSERT(real(u[1]) == 20.0);
   TEST_ASSERT(imag(u[2]) == 30.1);
   TEST_ASSERT(i2 == 13);
   TEST_ASSERT(u.capacity() == 3);
   #endif

}
#endif // if 0
 
TEST_BEGIN(GridArrayTest)
TEST_ADD(GridArrayTest, testConstructor)
TEST_ADD(GridArrayTest, testAllocate)
TEST_ADD(GridArrayTest, testSubscript)
//TEST_ADD(GridArrayTest, testCopyConstructor)
//TEST_ADD(GridArrayTest, testAssignment)
//TEST_ADD(GridArrayTest, testBaseClassReference)
//TEST_ADD(GridArrayTest, testSerializeFile1)
//TEST_ADD(GridArrayTest, testSerializeFile2)
TEST_END(GridArrayTest)

#endif
