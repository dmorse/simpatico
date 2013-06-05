#ifndef G_ARRAY_TEST_H
#define G_ARRAY_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <util/containers/GArray.h>
#include <util/containers/Array.h>

using namespace Util;

class GArrayTest : public UnitTest 
{

public:

   void setUp(){};
   void tearDown(){};
   void testReserve();
   void testConstructor();
   void testSubscript();
   void testDefaultReserve();
   void testResize();
   void testCopyConstructor();
   void testSerialize1File();

};


void GArrayTest::testConstructor()
{
   printMethod(TEST_FUNC);
   GArray<int> v;
   TEST_ASSERT(v.capacity() == 0);
} 

void GArrayTest::testReserve()
{
   printMethod(TEST_FUNC);
   GArray<int> v;
   v.reserve(3);
   TEST_ASSERT(v.capacity() == 3 );
   TEST_ASSERT(v.size() == 0 );
}

void GArrayTest::testSubscript()
{
   printMethod(TEST_FUNC);
   GArray<int> v;
   v.reserve(3);
   v.append(3);
   v.append(4);
   v.append(5);
   TEST_ASSERT(v.size() == 3);
   TEST_ASSERT(v[0] == 3);
   TEST_ASSERT(v[1] == 4);
   TEST_ASSERT(v[2] == 5);
} 

void GArrayTest::testDefaultReserve()
{
   printMethod(TEST_FUNC);
   GArray<int> v;
   v.append(3);
   v.append(4);
   v.append(5);
   TEST_ASSERT(v.size() == 3);
   TEST_ASSERT(v.capacity() == 64);
   TEST_ASSERT(v[0] == 3);
   TEST_ASSERT(v[1] == 4);
   TEST_ASSERT(v[2] == 5);
   v.append(6);
   v.append(7);
   TEST_ASSERT(v.size() == 5);
   TEST_ASSERT(v.capacity() == 64);
} 

void GArrayTest::testResize()
{
   printMethod(TEST_FUNC);
   GArray<int> v;
   v.reserve(3);
   v.append(3);
   v.append(4);
   v.append(5);
   TEST_ASSERT(v.size() == 3);
   TEST_ASSERT(v.capacity() == 3);
   TEST_ASSERT(v[0] == 3);
   TEST_ASSERT(v[1] == 4);
   TEST_ASSERT(v[2] == 5);
   v.append(6);
   v.append(7);
   TEST_ASSERT(v.size() == 5);
   TEST_ASSERT(v.capacity() == 6);
} 

void GArrayTest::testCopyConstructor()
{
   printMethod(TEST_FUNC);
   GArray<int> v;
   TEST_ASSERT(v.capacity() == 0 );

   v.reserve(3);
   TEST_ASSERT(v.capacity() == 3 );
   v.append(3);
   v.append(4);
   v.append(5);

   GArray<int> u(v);
   TEST_ASSERT(u.capacity() == 3 );
   TEST_ASSERT(u.size() == 3);
   TEST_ASSERT(v.size() == 3);
   TEST_ASSERT(v[0] == 3 );
   TEST_ASSERT(v[1] == 4 );
   TEST_ASSERT(v[2] == 5 );
   TEST_ASSERT(u[0] == 3 );
   TEST_ASSERT(u[1] == 4 );
   TEST_ASSERT(u[2] == 5 );
}


void GArrayTest::testSerialize1File()
{
   printMethod(TEST_FUNC);

   GArray<int> v;
   int i1 = 13;
   
   TEST_ASSERT(v.capacity() == 0 );

   v.reserve(4);
   TEST_ASSERT(v.capacity() == 4);
   v.append(3);
   v.append(4);
   v.append(5);

   TEST_ASSERT(v.size() == 3);
   TEST_ASSERT(v.capacity() == 4);
   TEST_ASSERT(v[0] == 3 );
   TEST_ASSERT(v[1] == 4 );
   TEST_ASSERT(v[2] == 5 );

   BinaryFileOArchive oArchive;
   openOutputFile("binary", oArchive.file());
   oArchive << v;
   oArchive << i1;
   oArchive.file().close();

   // Show that v is unchanged by packing
   TEST_ASSERT(v.size() == 3);
   TEST_ASSERT(v.capacity() == 4);
   TEST_ASSERT(v[0] == 3 );
   TEST_ASSERT(v[1] == 4 );
   TEST_ASSERT(v[2] == 5 );

   GArray<int> u;
   int i2;
   u.reserve(3);

   BinaryFileIArchive iArchive;
   openInputFile("binary", iArchive.file());
   iArchive >> u;
   iArchive >> i2;
   iArchive.file().close();

   TEST_ASSERT(u.size() == 3);
   TEST_ASSERT(u.capacity() == 4);
   TEST_ASSERT(u[0] == 3);
   TEST_ASSERT(u[1] == 4);
   TEST_ASSERT(u[2] == 5);
   TEST_ASSERT(i2 == 13);

}

TEST_BEGIN(GArrayTest)
TEST_ADD(GArrayTest, testReserve)
TEST_ADD(GArrayTest, testConstructor)
TEST_ADD(GArrayTest, testSubscript)
TEST_ADD(GArrayTest, testDefaultReserve)
TEST_ADD(GArrayTest, testResize)
TEST_ADD(GArrayTest, testCopyConstructor)
TEST_ADD(GArrayTest, testSerialize1File)
TEST_END(GArrayTest)

#endif
