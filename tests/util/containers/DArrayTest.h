#ifndef DARRAY_TEST_H
#define DARRAY_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <util/containers/DArray.h>
#include <util/archives/MemoryOArchive.h>
#include <util/archives/MemoryIArchive.h>
#include <util/archives/MemoryCounter.h>

using namespace Util;

class DArrayTest : public UnitTest 
{
private:

   const static int capacity = 3;

   typedef std::complex<double> Data;

   DArray<Data>  v;

public:

   void setUp() {}
   void tearDown() {}
   void testConstructor();
   void testAllocate();
   void testSubscript();
   void testIterator();
   void testCopyConstructor();
   void testAssignment();
   void testBaseClassReference();
   void testPack1();
   void testPack2();
   void testSerialize1();
   void testSerialize2();

};


void DArrayTest::testConstructor()
{
   printMethod(TEST_FUNC);
   TEST_ASSERT(v.capacity() == 0 );
   TEST_ASSERT(!v.isAllocated() );
} 

void DArrayTest::testAllocate()
{
   printMethod(TEST_FUNC);
   v.allocate(capacity);
   TEST_ASSERT(v.capacity() == capacity );
   TEST_ASSERT(v.isAllocated());
} 

void DArrayTest::testSubscript()
{
   printMethod(TEST_FUNC);
   v.allocate(capacity);
   for (int i=0; i < capacity; i++ ) {
      real(v[i]) = (i+1)*10 ;
      imag(v[i]) = (i+1)*10 + 0.1;
   }

   TEST_ASSERT(real(v[0]) == 10 );
   TEST_ASSERT(imag(v[1]) == 20.1 );
   TEST_ASSERT(real(v[2]) == 30 );
} 

void DArrayTest::testIterator()
{
   printMethod(TEST_FUNC);
   v.allocate(capacity);
   for (int i=0; i < capacity; i++ ) {
      real(v[i]) = (i+1)*10 ;
      imag(v[i]) = (i+1)*10 + 0.1;
   }

   ArrayIterator<Data> it;
   v.begin(it);
   TEST_ASSERT(imag(*it) == 10.1 );
   TEST_ASSERT(!it.isEnd());
   TEST_ASSERT(it.notEnd());
   ++it;
   TEST_ASSERT(real(*it) == 20 );
   TEST_ASSERT(!it.isEnd());
   TEST_ASSERT(it.notEnd());
   ++it;
   TEST_ASSERT(imag(*it) == 30.1 );
   ++it;
   TEST_ASSERT(it.isEnd());
   TEST_ASSERT(!it.notEnd());
} 

void DArrayTest::testCopyConstructor()
{
   printMethod(TEST_FUNC);
   TEST_ASSERT(v.capacity() == 0 );
   TEST_ASSERT(!v.isAllocated() );

   v.allocate(capacity);
   TEST_ASSERT(v.capacity() == capacity );
   TEST_ASSERT(v.isAllocated() );
   for (int i=0; i < capacity; i++ ) {
      real(v[i]) = (i+1)*10 ;
      imag(v[i]) = (i+1)*10 + 0.1;
   }

   DArray<Data> u(v);
   TEST_ASSERT(u.capacity() == capacity );
   TEST_ASSERT(u.isAllocated() );
   TEST_ASSERT(real(v[0]) == 10 );
   TEST_ASSERT(imag(v[1]) == 20.1 );
   TEST_ASSERT(real(v[2]) == 30 );
   TEST_ASSERT(real(u[0]) == 10 );
   TEST_ASSERT(imag(u[1]) == 20.1 );
   TEST_ASSERT(real(u[2]) == 30 );
} 

void DArrayTest::testAssignment()
{
   printMethod(TEST_FUNC);

   v.allocate(capacity);
   TEST_ASSERT(v.capacity() == 3 );
   TEST_ASSERT(v.isAllocated() );

   DArray<Data> u;
   u.allocate(3);
   TEST_ASSERT(u.capacity() == 3 );
   TEST_ASSERT(u.isAllocated() );

   for (int i=0; i < capacity; i++ ) {
      real(v[i]) = (i+1)*10 ;
      imag(v[i]) = (i+1)*10 + 0.1;
   }

   u  = v;

   TEST_ASSERT(u.capacity() == 3 );
   TEST_ASSERT(u.isAllocated() );
   TEST_ASSERT(real(v[0]) == 10 );
   TEST_ASSERT(imag(v[1]) == 20.1 );
   TEST_ASSERT(real(v[2]) == 30 );
   TEST_ASSERT(real(u[0]) == 10 );
   TEST_ASSERT(imag(u[1]) == 20.1 );
   TEST_ASSERT(real(u[2]) == 30 );
} 

void DArrayTest::testBaseClassReference()
{
   printMethod(TEST_FUNC);
   v.allocate(3);
   for (int i=0; i < capacity; i++ ) {
      real(v[i]) = (i+1)*10 ;
      imag(v[i]) = (i+1)*10 + 0.1;
   }
   
   Array<Data>& u = v;
   TEST_ASSERT(real(u[0]) == 10 );
   TEST_ASSERT(imag(u[1]) == 20.1 );
   TEST_ASSERT(real(u[2]) == 30 );
}

#if 0
void DArrayTest::testPack1()
{
   printMethod(TEST_FUNC);
   v.allocate(3);
   for (int i=0; i < capacity; i++ ) {
      real(v[i]) = (i+1)*10 ;
      imag(v[i]) = (i+1)*10 + 0.1;
   }
   int size = v.packedSize();
   
   char* begin;
   char* current;
   char* end;
   begin = new char[size + 4]; // make buffer too large
   end = begin + size + 4;
   current = begin;
   
   v.pack(current, end);
   TEST_ASSERT(current == begin + size);
   TEST_ASSERT(end == current + 4);

   // Show that v is unchanged by packing
   TEST_ASSERT(imag(v[0])==10.1);
   TEST_ASSERT(real(v[1])==20.0);
   TEST_ASSERT(imag(v[2])==30.1);
   TEST_ASSERT(v.capacity() == 3);

   DArray<Data> u;
   u.allocate(3);
   current = begin;
   u.unpack(current, end);
   TEST_ASSERT(current == begin + size);
   TEST_ASSERT(end == current + 4);
   TEST_ASSERT(imag(u[0])==10.1);
   TEST_ASSERT(real(u[1])==20.0);
   TEST_ASSERT(imag(u[2])==30.1);
   TEST_ASSERT(u.capacity() == 3);

}
 
void DArrayTest::testPack2()
{
   printMethod(TEST_FUNC);
   v.allocate(3);
   for (int i=0; i < capacity; i++ ) {
      real(v[i]) = (i+1)*10 ;
      imag(v[i]) = (i+1)*10 + 0.1;
   }
   int size = v.packedSize();
  
   PackedData buffer; 
   buffer.allocate(size + 4);

   buffer.beginPacking();
   v.pack(buffer);
   TEST_ASSERT(buffer.endAllocated() == buffer.cursor() + 4);

   // Show that v is unchanged by packing
   TEST_ASSERT(imag(v[0])==10.1);
   TEST_ASSERT(real(v[1])==20.0);
   TEST_ASSERT(imag(v[2])==30.1);
   TEST_ASSERT(v.capacity() == 3);

   DArray<Data> u;
   u.allocate(3);
   buffer.beginUnpacking();
   u.unpack(buffer);
   TEST_ASSERT(buffer.cursor() == buffer.begin() + size);
   TEST_ASSERT(buffer.endAllocated() == buffer.cursor() + 4);
   TEST_ASSERT(imag(u[0]) == 10.1);
   TEST_ASSERT(real(u[1]) == 20.0);
   TEST_ASSERT(imag(u[2]) == 30.1);
   TEST_ASSERT(u.capacity() == 3);

}
#endif
 
void DArrayTest::testSerialize1()
{
   printMethod(TEST_FUNC);
   v.allocate(3);
   for (int i=0; i < capacity; i++ ) {
      real(v[i]) = (i+1)*10 ;
      imag(v[i]) = (i+1)*10 + 0.1;
   }
   int size = memorySize(v);
  
   int i1 = 13;
   int i2;

   MemoryOArchive oArchive;
   oArchive.allocate(size + 12);

   oArchive << v;
   TEST_ASSERT(oArchive.cursor() == oArchive.begin() + size);
   oArchive << i1;

   // Show that v is unchanged by packing
   TEST_ASSERT(imag(v[0])==10.1);
   TEST_ASSERT(real(v[1])==20.0);
   TEST_ASSERT(imag(v[2])==30.1);
   TEST_ASSERT(v.capacity() == 3);

   DArray<Data> u;
   u.allocate(3);

   MemoryIArchive iArchive;
   iArchive = oArchive;
   TEST_ASSERT(iArchive.begin()  == oArchive.begin());
   TEST_ASSERT(iArchive.cursor() == iArchive.begin());

   // Load into u and i2
   iArchive >> u;
   TEST_ASSERT(iArchive.begin() == oArchive.begin());
   TEST_ASSERT(iArchive.end() == oArchive.cursor());
   TEST_ASSERT(iArchive.cursor() == iArchive.begin() + size);

   iArchive >> i2;
   TEST_ASSERT(iArchive.cursor() == iArchive.end());
   TEST_ASSERT(iArchive.begin() == oArchive.begin());
   TEST_ASSERT(iArchive.end() == oArchive.cursor());

   TEST_ASSERT(imag(u[0]) == 10.1);
   TEST_ASSERT(real(u[1]) == 20.0);
   TEST_ASSERT(imag(u[2]) == 30.1);
   TEST_ASSERT(i2 == 13);
   TEST_ASSERT(u.capacity() == 3);

   // Release
   iArchive.release();
   TEST_ASSERT(!iArchive.isAllocated());
   TEST_ASSERT(iArchive.begin() == 0);
   TEST_ASSERT(iArchive.cursor() == 0);
   TEST_ASSERT(iArchive.end() == 0);
   TEST_ASSERT(oArchive.cursor() == oArchive.begin() + size + sizeof(int));

   // Clear values of u and i2
   for (int i=0; i < capacity; i++ ) {
      real(u[i]) = 0.0;
      imag(u[i]) = 0.0;
   }
   i2 = 0;

   // Reload into u and i2
   iArchive = oArchive;
   iArchive >> u;
   TEST_ASSERT(iArchive.begin() == oArchive.begin());
   TEST_ASSERT(iArchive.end() == oArchive.cursor());
   TEST_ASSERT(iArchive.cursor() == iArchive.begin() + size);

   iArchive >> i2;
   TEST_ASSERT(iArchive.cursor() == iArchive.end());
   TEST_ASSERT(iArchive.begin() == oArchive.begin());
   TEST_ASSERT(iArchive.end() == oArchive.cursor());

   TEST_ASSERT(imag(u[0]) == 10.1);
   TEST_ASSERT(real(u[1]) == 20.0);
   TEST_ASSERT(imag(u[2]) == 30.1);
   TEST_ASSERT(i2 == 13);
   TEST_ASSERT(u.capacity() == 3);

}

void DArrayTest::testSerialize2()
{
   printMethod(TEST_FUNC);
   v.allocate(capacity);
   for (int i=0; i < capacity; i++ ) {
      real(v[i]) = (i+1)*10 ;
      imag(v[i]) = (i+1)*10 + 0.1;
   }
   int size = memorySize(v);
  
   MemoryOArchive oArchive;
   oArchive.allocate(size);

   oArchive << v;
   TEST_ASSERT(oArchive.cursor() == oArchive.begin() + size);

   // Show that v is unchanged by packing
   TEST_ASSERT(imag(v[0])==10.1);
   TEST_ASSERT(real(v[1])==20.0);
   TEST_ASSERT(imag(v[2])==30.1);
   TEST_ASSERT(v.capacity() == capacity);

   DArray<Data> u;

   // Note: We do not allocate DArray<Data> u in this test.
   // This is the main difference from testSerialize1()

   MemoryIArchive iArchive;

   iArchive = oArchive;

   TEST_ASSERT(iArchive.begin()  == oArchive.begin());
   TEST_ASSERT(iArchive.cursor() == iArchive.begin());

   iArchive >> u;

   TEST_ASSERT(iArchive.cursor() == iArchive.begin() + size);
   TEST_ASSERT(imag(u[0]) == 10.1);
   TEST_ASSERT(real(u[1]) == 20.0);
   TEST_ASSERT(imag(u[2]) == 30.1);
   TEST_ASSERT(u.capacity() == 3);

}

TEST_BEGIN(DArrayTest)
TEST_ADD(DArrayTest, testConstructor)
TEST_ADD(DArrayTest, testAllocate)
TEST_ADD(DArrayTest, testSubscript)
TEST_ADD(DArrayTest, testIterator)
TEST_ADD(DArrayTest, testCopyConstructor)
TEST_ADD(DArrayTest, testAssignment)
TEST_ADD(DArrayTest, testBaseClassReference)
//TEST_ADD(DArrayTest, testPack1)
//TEST_ADD(DArrayTest, testPack2)
TEST_ADD(DArrayTest, testSerialize1)
TEST_ADD(DArrayTest, testSerialize2)
TEST_END(DArrayTest)

#endif
