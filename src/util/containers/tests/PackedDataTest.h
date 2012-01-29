#ifndef PACKED_DATA_TEST_H
#define PACKED_DATA_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <util/containers/PackedData.h>

using namespace Util;

class PackedDataTest : public UnitTest 
{
private:

   const static int capacity = 64;

   PackedData  v;

public:

   void setUp() {}
   void tearDown() {}
   void testConstructor();
   void testAllocate();
   void testPack();
   void testPackArray();

};


void PackedDataTest::testConstructor()
{
   printMethod(TEST_FUNC);
   TEST_ASSERT(v.capacity() == 0 );
   TEST_ASSERT(!v.isAllocated() );
} 

void PackedDataTest::testAllocate()
{
   printMethod(TEST_FUNC);
   v.allocate(capacity);
   TEST_ASSERT(v.capacity() == capacity );
   TEST_ASSERT(v.isAllocated() );
} 

void PackedDataTest::testPack()
{
   printMethod(TEST_FUNC);
   v.allocate(capacity);
   TEST_ASSERT(v.endAllocated() == v.begin() + v.capacity());

   int offset;   
   int i1, i2;
   double d1, d2;
   std::complex<double> c1, c2;

   i1 = 3;
   d1 = 45.0;
   c1 = std::complex<double>(3.0, 4.0);
   
   v.beginPacking();
   TEST_ASSERT(v.cursor() == v.begin());
   v.pack(i1);
   offset = sizeof(int);
   TEST_ASSERT(v.cursor() == v.begin() + offset);
   v.pack(d1);
   offset += sizeof(double);
   TEST_ASSERT(v.cursor() == v.begin() + offset);
   v.pack(c1);
   offset += sizeof(std::complex<double>);
   TEST_ASSERT(v.cursor() == v.begin() + offset);

   v.beginUnpacking();
   TEST_ASSERT(v.endPacked() == v.begin() + offset);
   TEST_ASSERT(v.cursor() == v.begin());

   v.unpack(i2);
   offset = sizeof(int);
   TEST_ASSERT(v.cursor() == v.begin() + offset);
   v.unpack(d2);
   offset += sizeof(double);
   TEST_ASSERT(v.cursor() == v.begin() + offset);
   v.unpack(c2);
   offset += sizeof(std::complex<double>);
   TEST_ASSERT(v.cursor() == v.begin() + offset);

   TEST_ASSERT(i1 == i2);
   TEST_ASSERT(d1 == d2);
   TEST_ASSERT(c1 == c2);

}
 
void PackedDataTest::testPackArray()
{
   printMethod(TEST_FUNC);
   v.allocate(capacity);
   TEST_ASSERT(v.endAllocated() == v.begin() + v.capacity());
   TEST_ASSERT(v.endPacked() == v.begin());

   int offset;   
   int i1, i2;
   double d1, d2;
   std::complex<double> c1, c2;
   double a1[4];
   double a2[4];

   i1 = 3;
   d1 = 45.0;
   c1 = std::complex<double>(3.0, 4.0);
   a1[0] = 9.0;
   a1[1] = 8.0;
   a1[2] = 7.0;
   a1[3] = 6.0;
   
   v.beginPacking();
   TEST_ASSERT(v.cursor() == v.begin());
   TEST_ASSERT(v.endPacked() == v.begin());
   v.pack(i1);
   offset = sizeof(int);
   TEST_ASSERT(v.cursor() == v.begin() + offset);
   v.pack(d1);
   offset += sizeof(double);
   TEST_ASSERT(v.cursor() == v.begin() + offset);
   v.pack(a1, 4);
   offset += 4*sizeof(double);
   TEST_ASSERT(v.cursor() == v.begin() + offset);
   v.pack(c1);
   offset += sizeof(std::complex<double>);
   TEST_ASSERT(v.cursor() == v.begin() + offset);
   TEST_ASSERT(v.endPacked() == v.begin());

   v.beginUnpacking();
   TEST_ASSERT(v.cursor() == v.begin());
   TEST_ASSERT(v.endPacked() == v.begin() + offset);

   v.unpack(i2);
   offset = sizeof(int);
   TEST_ASSERT(v.cursor() == v.begin() + offset);
   v.unpack(d2);
   offset += sizeof(double);
   TEST_ASSERT(v.cursor() == v.begin() + offset);
   v.unpack(a2, 4);
   offset += 4*sizeof(double);
   TEST_ASSERT(v.cursor() == v.begin() + offset);
   v.unpack(c2);
   offset += sizeof(std::complex<double>);
   TEST_ASSERT(v.cursor() == v.begin() + offset);

   TEST_ASSERT(i1 == i2);
   TEST_ASSERT(d1 == d2);
   TEST_ASSERT(c1 == c2);
   for (int j = 0; j < 4; ++j) {
      TEST_ASSERT(a1[j] == a2[j]);
   }

}

TEST_BEGIN(PackedDataTest)
TEST_ADD(PackedDataTest, testConstructor)
TEST_ADD(PackedDataTest, testAllocate)
TEST_ADD(PackedDataTest, testPack)
TEST_ADD(PackedDataTest, testPackArray)
TEST_END(PackedDataTest)

#endif
