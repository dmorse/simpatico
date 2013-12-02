#ifndef DDMD_ATOM_MAP_TEST_H
#define DDMD_ATOM_MAP_TEST_H

#include <ddMd/storage/AtomMap.h>
#include <ddMd/chemistry/Atom.h>
#include <ddMd/chemistry/Group.h>
#include <ddMd/chemistry/AtomArray.h>
#include <util/boundary/OrthorhombicBoundary.h>
#include <util/space/Vector.h>
#include <util/containers/ArraySet.h>

#ifdef UTIL_MPI
#ifndef TEST_MPI
#define TEST_MPI
#endif
#endif

#include <test/ParamFileTest.h>
#include <test/UnitTestRunner.h>

using namespace Util;
using namespace DdMd;

class AtomMapTest : public ParamFileTest
{

private:

   AtomArray array_;
   AtomMap map_;
   int n_;

public:

   virtual void setUp()
   { 
      int n_ = 20;
      array_.allocate(n_);
      map_.allocate(n_);
      int bit = 0;
      for (int i=0; i < n_; ++i) {
         array_[i].setId(i);
         if (bit) {
            array_[i].setIsGhost(true);
            bit = 0;
         } else {
            array_[i].setIsGhost(false);
            bit = 1;
         }
      }
   }

   void testAdd();
   void testAddRemove();
   void testFindLocal();
   void testFindLocalGhost();
   void testClearGhosts();
   void testFindGroupGhostAtoms1();
   void testFindGroupGhostAtoms2();
   void testFindGroupGhostAtoms3();

};

inline void AtomMapTest::testAdd()
{
   printMethod(TEST_FUNC);

   array_[13].setId(11);
   array_[8].setId(6);

   map_.addLocal(&array_[6]);
   //map_.addLocal(&array_[8]);
   map_.addLocal(&array_[12]);
   map_.addGhost(&array_[9]);
   map_.addGhost(&array_[11]);
   TEST_ASSERT(map_.nLocal() == 2);
   TEST_ASSERT(map_.nGhostDistinct() == 2);
   TEST_ASSERT(map_.isValid());
   map_.addLocal(&array_[14]);
   map_.addGhost(&array_[13]);
   TEST_ASSERT(map_.nLocal() == 3);
   TEST_ASSERT(map_.nGhostDistinct() == 2);
   TEST_ASSERT(map_.isValid());
}

inline void AtomMapTest::testAddRemove()
{
   printMethod(TEST_FUNC);

   array_[13].setId(11);
   array_[8].setId(6);

   map_.addLocal(&array_[6]);
   //map_.addLocal(&array_[8]);
   map_.addLocal(&array_[12]);
   map_.addGhost(&array_[9]);
   map_.addLocal(&array_[14]);
   map_.addGhost(&array_[11]);
   map_.addGhost(&array_[13]);
   TEST_ASSERT(map_.nLocal() == 3);
   TEST_ASSERT(map_.nGhostDistinct() == 2);
   TEST_ASSERT(map_.isValid());

   map_.removeLocal(&array_[6]);
   map_.removeGhost(&array_[13]);
   map_.removeLocal(&array_[12]);
   map_.removeGhost(&array_[11]);
   TEST_ASSERT(map_.nLocal() == 1);
   TEST_ASSERT(map_.nGhostDistinct() == 1);
   TEST_ASSERT(map_.isValid());

}

inline void AtomMapTest::testFindLocal()
{
   printMethod(TEST_FUNC);

   array_[13].setId(11);
   array_[8].setId(6);

   map_.addLocal(&array_[6]);
   //map_.addLocal(&array_[8]);
   map_.addLocal(&array_[12]);
   map_.addGhost(&array_[9]);
   map_.addLocal(&array_[14]);
   map_.addGhost(&array_[11]);
   map_.addGhost(&array_[13]);
   TEST_ASSERT(map_.nLocal() == 3);
   TEST_ASSERT(map_.nGhostDistinct() == 2);
   TEST_ASSERT(map_.isValid());

   TEST_ASSERT( array_[14].id() == 14);
   TEST_ASSERT( map_.find(14) == &array_[14] );
   TEST_ASSERT( map_.find(12) == &array_[12] );
   TEST_ASSERT( map_.find(6) == &array_[6] );

   map_.removeLocal(&array_[6]);
   map_.removeGhost(&array_[13]);
   TEST_ASSERT(map_.isValid());
   map_.removeLocal(&array_[12]);
   map_.removeGhost(&array_[11]);
   TEST_ASSERT(map_.isValid());
   TEST_ASSERT(map_.nLocal() == 1);
   TEST_ASSERT(map_.nGhostDistinct() == 1);
   TEST_ASSERT( map_.find(6) == 0 );
   TEST_ASSERT( map_.find(14) == &array_[14] );
   TEST_ASSERT( map_.find(12) == 0 );

}

inline void AtomMapTest::testFindLocalGhost()
{
   printMethod(TEST_FUNC);

   array_[13].setId(11);
   array_[15].setId(11);
   array_[8].setId(6);
   array_[7].setId(6);

   map_.addLocal(&array_[6]);
   map_.addLocal(&array_[12]);
   map_.addGhost(&array_[9]);
   map_.addLocal(&array_[14]);
   map_.addGhost(&array_[11]);
   TEST_ASSERT(map_.nGhost() == 2);
   TEST_ASSERT(map_.nGhostDistinct() == 2);
   map_.addGhost(&array_[13]);
   TEST_ASSERT(map_.nGhostDistinct() == 2);
   TEST_ASSERT(map_.nGhost() == 3);
   map_.addGhost(&array_[15]);
   TEST_ASSERT(map_.nGhostDistinct() == 2);
   TEST_ASSERT(map_.nGhost() == 4);
   TEST_ASSERT(map_.nLocal() == 3);
   TEST_ASSERT(map_.isValid());

   TEST_ASSERT( array_[14].id() == 14);
   TEST_ASSERT( map_.find(14) == &array_[14] );
   TEST_ASSERT( map_.find(12) == &array_[12] );
   TEST_ASSERT( map_.find(6) == &array_[6] );
   TEST_ASSERT( map_.find(9) == &array_[9] );
   TEST_ASSERT( map_.find(11) == &array_[11] );
   TEST_ASSERT( map_.find(15) == 0 );
   TEST_ASSERT( map_.find(16) == 0 );
   TEST_ASSERT( map_.find(17) == 0 );

   map_.addGhost(&array_[7]); 
   TEST_ASSERT(map_.nGhost() == 5);
   TEST_ASSERT(map_.nGhostDistinct() == 2);
   TEST_ASSERT(map_.nLocal() == 3);
   TEST_ASSERT(map_.find(6) == &array_[6]);

   map_.removeLocal(&array_[6]);
   map_.removeGhost(&array_[13]);
   TEST_ASSERT(map_.nGhost() == 4);
   TEST_ASSERT(map_.nGhostDistinct() == 3);
   TEST_ASSERT(map_.find(6) == &array_[7]);
   TEST_ASSERT(map_.find(9) == &array_[9]);
   TEST_ASSERT(map_.find(11) == &array_[11]);
   map_.removeLocal(&array_[12]);
   map_.removeGhost(&array_[11]);
   TEST_ASSERT(map_.nGhost() == 3);
   TEST_ASSERT(map_.nGhostDistinct() == 3);
   TEST_ASSERT(map_.find(11) == &array_[15]);
   TEST_ASSERT(map_.nGhost() == 3);
   TEST_ASSERT(map_.nGhostDistinct() == 3);
   TEST_ASSERT(map_.nLocal() == 1);
   TEST_ASSERT(map_.isValid());
   TEST_ASSERT(map_.find(6) == &array_[7]);
   TEST_ASSERT(map_.find(14) == &array_[14] );
   TEST_ASSERT(map_.find(12) == 0);
   TEST_ASSERT(map_.find(3) == 0);
   TEST_ASSERT(map_.find(5) == 0);
   TEST_ASSERT(map_.find(7) == 0);
   TEST_ASSERT(map_.find(9) == &array_[9] );
   TEST_ASSERT(map_.find(13) == 0);

}

inline void AtomMapTest::testClearGhosts()
{
   printMethod(TEST_FUNC);
   ArraySet<Atom> ghostSet;
   ghostSet.allocate(array_);

   array_[13].setId(11);
   array_[15].setId(11);
   array_[8].setId(6);
   array_[7].setId(6);

   map_.addLocal(&array_[6]);
   map_.addLocal(&array_[12]);
   map_.addGhost(&array_[9]); ghostSet.append(array_[9]);
   map_.addLocal(&array_[14]);
   map_.addGhost(&array_[11]); ghostSet.append(array_[11]);
   TEST_ASSERT(map_.nGhost() == 2);
   TEST_ASSERT(map_.nGhostDistinct() == 2);
   map_.addGhost(&array_[13]); ghostSet.append(array_[13]);
   TEST_ASSERT(map_.nGhost() == 3);
   TEST_ASSERT(map_.nGhostDistinct() == 2);
   map_.addGhost(&array_[15]); ghostSet.append(array_[15]);
   TEST_ASSERT(map_.nGhost() == 4);
   TEST_ASSERT(map_.nGhostDistinct() == 2);
   TEST_ASSERT(map_.nLocal() == 3);
   TEST_ASSERT(map_.isValid());

   
   TEST_ASSERT(map_.nGhost() == 4);
   TEST_ASSERT(map_.nGhostDistinct() == 2);
   map_.removeLocal(&array_[6]);
   TEST_ASSERT(map_.nLocal() == 2);
   map_.removeGhost(&array_[13]); ghostSet.remove(array_[13]);
   TEST_ASSERT(map_.nGhost() == 3);
   TEST_ASSERT(map_.nGhostDistinct() == 2);
   TEST_ASSERT(map_.isValid());

   map_.addGhost(&array_[1]); ghostSet.append(array_[1]);
   map_.addGhost(&array_[3]); ghostSet.append(array_[3]);
   TEST_ASSERT(map_.nGhost() == 5);
   TEST_ASSERT(map_.nGhostDistinct() == 4);

   map_.clearGhosts(ghostSet);
   TEST_ASSERT(map_.isValid());
}

inline void AtomMapTest::testFindGroupGhostAtoms1()
{
   printMethod(TEST_FUNC);

   Vector L(1.0, 1.0, 1.0);
   OrthorhombicBoundary boundary;
   boundary.setOrthorhombic(L);

   Group<2> bond;
   bond.setAtomId(0, 4);
   bond.setAtomId(1, 2);
   bond.setAtomPtr(0, &array_[4]);
   bond.setAtomPtr(1, &array_[2]);
   array_[4].position() = Vector(0.9, 0.4, 0.4);
   array_[2].position() = Vector(0.7, 0.5, 0.5);

   // Both atoms local, satisfy minimum image convention.
   map_.addLocal(&array_[2]);
   map_.addLocal(&array_[4]);
   TEST_ASSERT(map_.nLocal() == 2);
   TEST_ASSERT(map_.nGhost() == 0);
   TEST_ASSERT(map_.nGhostDistinct() == 0);
   TEST_ASSERT(map_.isValid());

   map_.findGroupGhostAtoms(bond, boundary);
}

inline void AtomMapTest::testFindGroupGhostAtoms2()
{
   printMethod(TEST_FUNC);

   Vector L(1.0, 1.0, 1.0);
   OrthorhombicBoundary boundary;
   boundary.setOrthorhombic(L);

   array_[3].setId(2);
   array_[3].position() = Vector(1.1, 0.5, 0.5);
   array_[4].position() = Vector(0.9, 0.4, 0.4);
   TEST_ASSERT(array_[4].id() == 4);
   TEST_ASSERT(array_[3].id() == 2);

   Group<2> bond;
   bond.setAtomId(0, 2); // ghost
   bond.setAtomId(1, 4); // local
   bond.setAtomPtr(1, &array_[4]);
   // bond.setAtomPtr(1, &array_[3]);

   // Both atoms local, satisfy minimum image convention.
   map_.addLocal(&array_[4]);
   map_.addGhost(&array_[3]);
   TEST_ASSERT(map_.nLocal() == 1);
   TEST_ASSERT(map_.nGhost() == 1);
   TEST_ASSERT(map_.nGhostDistinct() == 1);
   TEST_ASSERT(map_.isValid());

   map_.findGroupGhostAtoms(bond, boundary);
   TEST_ASSERT(bond.atomPtr(0) == &array_[3]);
   TEST_ASSERT(bond.atomPtr(1) == &array_[4]);
}

inline void AtomMapTest::testFindGroupGhostAtoms3()
{
   printMethod(TEST_FUNC);

   Vector L(1.0, 1.0, 1.0);
   OrthorhombicBoundary boundary;
   boundary.setOrthorhombic(L);

   array_[2].position() = Vector(0.1, 0.5, 0.5);
   array_[4].position() = Vector(0.9, 0.4, 0.4);
   array_[3].setId(2);
   array_[3].position() = Vector(1.1, 0.5, 0.5);  // ghost of 2
   array_[5].setId(4);
   array_[5].position() = Vector(-0.1, 0.4, 0.4); // ghost of 4
   TEST_ASSERT(array_[2].id() == 2);
   TEST_ASSERT(array_[3].id() == 2);
   TEST_ASSERT(array_[4].id() == 4);
   TEST_ASSERT(array_[5].id() == 4);

   Group<2> bond;
   bond.setAtomId(0, 2); // local and ghost
   bond.setAtomId(1, 4); // local
   bond.setAtomPtr(0, &array_[2]);
   bond.setAtomPtr(1, &array_[4]);

   // One atom local only, one local and ghost.
   map_.addLocal(&array_[2]);
   map_.addLocal(&array_[4]);
   map_.addGhost(&array_[3]);
   map_.addGhost(&array_[5]);
   TEST_ASSERT(map_.nLocal() == 2);
   TEST_ASSERT(map_.nGhost() == 2);
   TEST_ASSERT(map_.nGhostDistinct() == 0);
   TEST_ASSERT(map_.isValid());

   map_.findGroupGhostAtoms(bond, boundary);
   TEST_ASSERT(bond.atomPtr(0) == &array_[2]);
   TEST_ASSERT(bond.atomPtr(1) == &array_[5]);
}

TEST_BEGIN(AtomMapTest)
TEST_ADD(AtomMapTest, testAdd)
TEST_ADD(AtomMapTest, testAddRemove)
TEST_ADD(AtomMapTest, testFindLocal)
TEST_ADD(AtomMapTest, testFindLocalGhost)
TEST_ADD(AtomMapTest, testClearGhosts)
TEST_ADD(AtomMapTest, testFindGroupGhostAtoms1)
TEST_ADD(AtomMapTest, testFindGroupGhostAtoms2)
TEST_ADD(AtomMapTest, testFindGroupGhostAtoms3)
TEST_END(AtomMapTest)

#endif
