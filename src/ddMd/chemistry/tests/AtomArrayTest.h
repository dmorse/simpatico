#ifndef ATOM_ARRAY_TEST_H
#define ATOM_ARRAY_TEST_H

#include <ddMd/chemistry/AtomArray.h>
#include <ddMd/chemistry/Atom.h>
#include <util/space/Vector.h>

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

using namespace Util;
using namespace DdMd;

class AtomArrayTest : public UnitTest 
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


void AtomArrayTest::testConstructor()
{
   printMethod(TEST_FUNC);
   AtomArray v;
   TEST_ASSERT(v.capacity() == 0 );
   TEST_ASSERT(!v.isAllocated() );
} 

void AtomArrayTest::testAllocate()
{
   printMethod(TEST_FUNC);
   AtomArray v;
   v.allocate(3);
   TEST_ASSERT(v.capacity() == 3 );
   TEST_ASSERT(v.isAllocated() );
} 

void AtomArrayTest::testSubscript()
{
   printMethod(TEST_FUNC);
   AtomArray a;
   a.allocate(3);
   Vector r, v, f;

   r[0] = 3.0;
   r[1] = 4.0;
   r[2] = 5.0;

   v[0] = -1.0;
   v[1] =  2.5;
   v[2] = -0.5;

   f[0] =  2.0;
   f[1] =  0.2;
   f[2] = -8.5;

   TEST_ASSERT(!a[1].isGhost());

   a[1].setId(1);
   a[1].setTypeId(5);
   a[1].position() = r;
   a[1].velocity() = v;
   a[1].force()    = f;
   a[1].setIsGhost(true);

   TEST_ASSERT(a[1].id() == 1 );
   TEST_ASSERT(a[1].typeId() == 5 );
   TEST_ASSERT(a[1].position() == r );
   TEST_ASSERT(a[1].velocity() == v );
   TEST_ASSERT(a[1].force()    == f );
   TEST_ASSERT(a[1].isGhost());

   a[2].setId(1);
   a[2].setTypeId(5);
   a[2].position() = f;
   a[2].velocity() = r;
   a[2].force()    = v;
   TEST_ASSERT(a[2].position() == a[1].force() );
   TEST_ASSERT(a[2].velocity() == a[1].position() );
   TEST_ASSERT(a[2].force()    == a[1].velocity() );

} 

TEST_BEGIN(AtomArrayTest)
TEST_ADD(AtomArrayTest, testConstructor)
TEST_ADD(AtomArrayTest, testAllocate)
TEST_ADD(AtomArrayTest, testSubscript)
TEST_END(AtomArrayTest)

#endif
