#ifndef DDMD_SP_ATOM_STORAGE_TEST_H
#define DDMD_SP_ATOM_STORAGE_TEST_H

#include <ddMd/sp/storage/SpAtomStorage.h>
#include <ddMd/sp/chemistry/SpAtom.h>
#include <ddMd/sp/chemistry/SpGroup.h>

#include <test/ParamFileTest.h>
#include <test/UnitTestRunner.h>

using namespace Util;
using namespace DdMd;

class SpAtomStorageTest : public ParamFileTest
{

private:

   SpAtomStorage storage_;

public:

   SpAtomStorageTest() 
    : storage_()
   {}

   // virtual void setUp() { }

   void testAddAtoms();

};

inline void SpAtomStorageTest::testAddAtoms()
{
   printMethod(TEST_FUNC);

   storage_.allocate(100);

   SpAtom* atomPtr22;
   atomPtr22 = storage_.newPtr();
   atomPtr22->id = 22;
   storage_.add();
   TEST_ASSERT(atomPtr22 == storage_.ptr(22));
   TEST_ASSERT(0 == storage_.ptr(15));
   TEST_ASSERT(0 == storage_.ptr(13));
   TEST_ASSERT(storage_.size() == 1);

   SpAtom* atomPtr13;
   atomPtr13 = storage_.newPtr();
   atomPtr13->id = 13;
   storage_.add();
   TEST_ASSERT(atomPtr13 == storage_.ptr(13));
   TEST_ASSERT(atomPtr22 == storage_.ptr(22));
   TEST_ASSERT(0 == storage_.ptr(15));
   TEST_ASSERT(storage_.size() == 2);

   SpAtom* atomPtr15;
   atomPtr15 = storage_.newPtr();
   atomPtr15->id = 15;
   storage_.add();
   TEST_ASSERT(atomPtr15 == storage_.ptr(15));
   TEST_ASSERT(atomPtr13 == storage_.ptr(13));
   TEST_ASSERT(atomPtr22 == storage_.ptr(22));
   TEST_ASSERT(storage_.size() == 3);
}

TEST_BEGIN(SpAtomStorageTest)
TEST_ADD(SpAtomStorageTest, testAddAtoms)
TEST_END(SpAtomStorageTest)

#endif
