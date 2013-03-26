#ifndef DDMD_ATOM_TEST_H
#define DDMD_ATOM_TEST_H

#include <ddMd/chemistry/AtomArray.h>
#include <ddMd/chemistry/Atom.h>
#include <util/space/Vector.h>

#ifdef UTIL_MPI
#ifndef TEST_MPI
#define TEST_MPI
#endif
#endif

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

using namespace Util;
using namespace DdMd;

class AtomTest : public UnitTest 
{

public:

   void setUp() {}

   void tearDown() {}
  
   void testConstructor();

   void testAllocate();

   void testSubscript();

   void testAssignment();

};


void AtomTest::testConstructor()
{
   printMethod(TEST_FUNC);
   AtomArray v;
   TEST_ASSERT(v.capacity() == 0 );
   TEST_ASSERT(!v.isAllocated() );
} 

void AtomTest::testAllocate()
{
   printMethod(TEST_FUNC);
   AtomArray v;
   v.allocate(3);
   TEST_ASSERT(v.capacity() == 3 );
   TEST_ASSERT(v.isAllocated() );
} 

void AtomTest::testSubscript()
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

   TEST_ASSERT(!a[0].isGhost());
   TEST_ASSERT(!a[1].isGhost());
   TEST_ASSERT(!a[2].isGhost());
   TEST_ASSERT(a[0].id() == -1);
   TEST_ASSERT(a[1].id() == -1);
   TEST_ASSERT(a[2].id() == -1);

   a[1].setId(25);
   a[1].setTypeId(5);
   a[1].position() = r;
   a[1].velocity() = v;
   a[1].force()    = f;
   a[1].setIsGhost(true);
   a[1].mask().clear();
   a[1].mask().append(39);
   a[1].mask().append(37);
   a[1].plan().setFlags(23);

   TEST_ASSERT(a[1].id() == 25);
   TEST_ASSERT(a[1].typeId() == 5);
   TEST_ASSERT(a[1].position() == r);
   TEST_ASSERT(a[1].velocity() == v);
   TEST_ASSERT(a[1].force()    == f);
   TEST_ASSERT(a[1].isGhost());
   TEST_ASSERT(a[1].mask().isMasked(39));
   TEST_ASSERT(a[1].mask().isMasked(37));
   TEST_ASSERT(a[1].mask().size() == 2);
   TEST_ASSERT(!a[1].mask().isMasked(50));
   TEST_ASSERT(a[1].plan().flags() == 23);

   a[1].setIsGhost(false);
   TEST_ASSERT(a[1].typeId() == 5);
   TEST_ASSERT(a[1].position() == r);
   TEST_ASSERT(a[1].velocity() == v);
   TEST_ASSERT(a[1].force()    == f);
   TEST_ASSERT(!a[1].isGhost());
   TEST_ASSERT(a[1].mask().isMasked(39));
   TEST_ASSERT(a[1].mask().isMasked(37));
   TEST_ASSERT(!a[1].mask().isMasked(50));
   TEST_ASSERT(a[1].id() == 25);
   TEST_ASSERT(a[1].plan().flags() == 23);

   a[2].setId(35);
   a[2].setTypeId(8);
   a[2].position() = f;
   a[2].force()    = v;
   a[2].velocity() = r;
   a[2].setIsGhost(false);
   a[2].mask().clear();
   a[2].mask().append(59);
   a[2].mask().append(61);
   a[2].plan().setFlags(17);

   TEST_ASSERT(a[2].typeId() == 8);
   TEST_ASSERT(a[2].position() == a[1].force() );
   TEST_ASSERT(a[2].force()    == a[1].velocity() );
   TEST_ASSERT(a[2].velocity() == a[1].position() );
   TEST_ASSERT(!(a[2].isGhost()));
   TEST_ASSERT(a[2].mask().isMasked(59));
   TEST_ASSERT(a[2].mask().isMasked(61));
   TEST_ASSERT(!a[2].mask().isMasked(39));
   TEST_ASSERT(a[2].plan().flags() == 17);
   TEST_ASSERT(a[2].id() == 35);

   a[2].setIsGhost(true);
   TEST_ASSERT(a[2].position() == a[1].force() );
   TEST_ASSERT(a[2].force()    == a[1].velocity() );
   TEST_ASSERT(a[2].velocity() == a[1].position() );
   TEST_ASSERT(!a[1].isGhost());
   TEST_ASSERT(a[2].isGhost());
   TEST_ASSERT(a[2].velocity() == a[1].position() );
   TEST_ASSERT(a[2].mask().isMasked(59));
   TEST_ASSERT(a[2].mask().isMasked(61));
   TEST_ASSERT(!a[2].mask().isMasked(39));
   TEST_ASSERT(a[2].plan().flags() == 17);
   TEST_ASSERT(a[2].id() == 35);

} 

void AtomTest::testAssignment()
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

   TEST_ASSERT(!a[0].isGhost());
   TEST_ASSERT(!a[1].isGhost());
   TEST_ASSERT(!a[2].isGhost());
   TEST_ASSERT(a[0].id() == -1);
   TEST_ASSERT(a[1].id() == -1);
   TEST_ASSERT(a[2].id() == -1);
   TEST_ASSERT(a[0].mask().size() == 0);
   TEST_ASSERT(a[1].mask().size() == 0);
   TEST_ASSERT(a[2].mask().size() == 0);

   a[1].setId(25);
   a[1].setTypeId(5);
   a[1].position() = r;
   a[1].velocity() = v;
   a[1].force()    = f;
   a[1].setIsGhost(true);
   a[1].mask().clear();
   a[1].mask().append(39);
   a[1].mask().append(37);
   a[1].plan().setFlags(23);

   TEST_ASSERT(a[1].id() == 25);
   TEST_ASSERT(a[1].typeId() == 5);
   TEST_ASSERT(a[1].position() == r);
   TEST_ASSERT(a[1].velocity() == v);
   TEST_ASSERT(a[1].force()    == f);
   TEST_ASSERT(a[1].isGhost());
   TEST_ASSERT(a[1].mask().isMasked(39));
   TEST_ASSERT(a[1].mask().isMasked(37));
   TEST_ASSERT(a[1].mask().size() == 2);
   TEST_ASSERT(!a[1].mask().isMasked(50));
   TEST_ASSERT(a[1].plan().flags() == 23);

   a[1].setIsGhost(false);
   TEST_ASSERT(a[1].typeId() == 5);
   TEST_ASSERT(a[1].position() == r);
   TEST_ASSERT(a[1].force()    == f);
   TEST_ASSERT(!a[1].isGhost());
   TEST_ASSERT(a[1].velocity() == v);
   TEST_ASSERT(a[1].mask().isMasked(39));
   TEST_ASSERT(a[1].mask().isMasked(37));
   TEST_ASSERT(!a[1].mask().isMasked(50));
   TEST_ASSERT(a[1].plan().flags() == 23);
   TEST_ASSERT(a[1].id() == 25);

   a[2] = a[1];
   TEST_ASSERT(a[2].typeId() == 5);
   TEST_ASSERT(a[2].position() == r);
   TEST_ASSERT(a[2].force()    == f);
   TEST_ASSERT(!a[2].isGhost());
   TEST_ASSERT(a[2].velocity() == v);
   TEST_ASSERT(a[2].mask().isMasked(39));
   TEST_ASSERT(a[2].mask().isMasked(37));
   TEST_ASSERT(!a[2].mask().isMasked(50));
   TEST_ASSERT(a[2].plan().flags() == 23);
   TEST_ASSERT(a[2].id() == 25);

   a[2].setId(28);
   a[2].setIsGhost(true);
   TEST_ASSERT(a[2].id() == 28);
   TEST_ASSERT(a[2].isGhost());
   TEST_ASSERT(a[2].typeId() == 5);
   TEST_ASSERT(a[2].position() == r);
   TEST_ASSERT(a[2].force()    == f);
   TEST_ASSERT(a[2].velocity() == v);
   TEST_ASSERT(a[2].mask().isMasked(39));
   TEST_ASSERT(a[2].mask().isMasked(37));
   TEST_ASSERT(!a[2].mask().isMasked(50));
   TEST_ASSERT(a[2].plan().flags() == 23);

} 
TEST_BEGIN(AtomTest)
TEST_ADD(AtomTest, testConstructor)
TEST_ADD(AtomTest, testAllocate)
TEST_ADD(AtomTest, testSubscript)
TEST_ADD(AtomTest, testAssignment)
TEST_END(AtomTest)

#endif
