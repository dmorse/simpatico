#ifndef SPAN_ATOM_STORAGE_TEST_H
#define SPAN_ATOM_STORAGE_TEST_H

#include <spAn/storage/AtomStorage.h>
#include <spAn/chemistry/Atom.h>
#include <spAn/chemistry/Group.h>

#include <test/ParamFileTest.h>
#include <test/UnitTestRunner.h>

using namespace Util;
using namespace SpAn;

class AtomStorageTest : public ParamFileTest
{

private:

   AtomStorage storage_;

public:

   AtomStorageTest() 
    : storage_()
   {}

   // virtual void setUp() { }

   void testAddAtoms();

};

inline void AtomStorageTest::testAddAtoms()
{
   printMethod(TEST_FUNC);

   storage_.allocate(100);

   Atom* atomPtr22;
   atomPtr22 = storage_.newPtr();
   atomPtr22->id = 22;
   storage_.add();
   TEST_ASSERT(atomPtr22 == storage_.ptr(22));
   TEST_ASSERT(0 == storage_.ptr(15));
   TEST_ASSERT(0 == storage_.ptr(13));
   TEST_ASSERT(storage_.size() == 1);

   Atom* atomPtr13;
   atomPtr13 = storage_.newPtr();
   atomPtr13->id = 13;
   storage_.add();
   TEST_ASSERT(atomPtr13 == storage_.ptr(13));
   TEST_ASSERT(atomPtr22 == storage_.ptr(22));
   TEST_ASSERT(0 == storage_.ptr(15));
   TEST_ASSERT(storage_.size() == 2);

   Atom* atomPtr15;
   atomPtr15 = storage_.newPtr();
   atomPtr15->id = 15;
   storage_.add();
   TEST_ASSERT(atomPtr15 == storage_.ptr(15));
   TEST_ASSERT(atomPtr13 == storage_.ptr(13));
   TEST_ASSERT(atomPtr22 == storage_.ptr(22));
   TEST_ASSERT(storage_.size() == 3);
}

TEST_BEGIN(AtomStorageTest)
TEST_ADD(AtomStorageTest, testAddAtoms)
TEST_END(AtomStorageTest)

#endif
