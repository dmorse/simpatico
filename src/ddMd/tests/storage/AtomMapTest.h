#ifndef DDMD_ATOM_MAP_TEST_H
#define DDMD_ATOM_MAP_TEST_H

#include <ddMd/storage/AtomMap.h>
#include <ddMd/chemistry/Atom.h>
#include <ddMd/chemistry/AtomArray.h>
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
      n_ = 20;
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

TEST_BEGIN(AtomMapTest)
TEST_ADD(AtomMapTest, testAdd)
TEST_ADD(AtomMapTest, testAddRemove)
TEST_ADD(AtomMapTest, testFindLocal)
TEST_ADD(AtomMapTest, testFindLocalGhost)
TEST_ADD(AtomMapTest, testClearGhosts)
TEST_END(AtomMapTest)

#endif
