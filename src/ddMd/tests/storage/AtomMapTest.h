#ifndef DDMD_ATOM_MAP_TEST_H
#define DDMD_ATOM_MAP_TEST_H

#include <ddMd/storage/AtomMap.h>
#include <ddMd/chemistry/Atom.h>
#include <ddMd/chemistry/AtomArray.h>

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
      for (int i=0; i < n_; ++i) {
         array_[i].setId(i);
      }
   }

   void testAdd();
   void testAddRemove();

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
   UTIL_ASSERT(map_.nLocal() == 2);
   UTIL_ASSERT(map_.nGhost() == 2);
   UTIL_ASSERT(map_.isValid());
   map_.addLocal(&array_[14]);
   map_.addGhost(&array_[13]);
   UTIL_ASSERT(map_.nLocal() == 3);
   UTIL_ASSERT(map_.nGhost() == 2);
   UTIL_ASSERT(map_.isValid());
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
   UTIL_ASSERT(map_.nLocal() == 3);
   UTIL_ASSERT(map_.nGhost() == 2);
   UTIL_ASSERT(map_.isValid());

   map_.removeLocal(&array_[6]);
   map_.removeGhost(&array_[13]);
   map_.removeLocal(&array_[12]);
   map_.removeGhost(&array_[11]);
   UTIL_ASSERT(map_.nLocal() == 1);
   UTIL_ASSERT(map_.nGhost() == 1);
   UTIL_ASSERT(map_.isValid());

}
TEST_BEGIN(AtomMapTest)
TEST_ADD(AtomMapTest, testAdd)
TEST_ADD(AtomMapTest, testAddRemove)
TEST_END(AtomMapTest)

#endif
