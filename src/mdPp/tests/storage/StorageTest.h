#ifndef MDPP_STORAGE_TEST_H
#define MDPP_STORAGE_TEST_H

#include <mdPp/storage/Storage.h>
#include <mdPp/chemistry/Atom.h>
#include <mdPp/chemistry/Group.h>

#include <test/ParamFileTest.h>
#include <test/UnitTestRunner.h>

using namespace Util;
using namespace MdPp;

class StorageTest : public ParamFileTest
{

private:

   Storage storage_;

public:

   StorageTest() 
    : storage_()
   {}

   virtual void setUp()
   { 
      storage_.readParam("in/Storage"); 
   }

   void testReadParam();

   void testAddAtoms();

};

inline void StorageTest::testReadParam()
{  
   printMethod(TEST_FUNC); 
   if (verbose() > 0) {
      storage_.writeParam(std::cout);
   }
   TEST_ASSERT(storage_.nSpecies() == 0);
}

inline void StorageTest::testAddAtoms()
{
   printMethod(TEST_FUNC);

   Atom* atomPtr22;
   atomPtr22 = storage_.newAtomPtr();
   atomPtr22->id = 22;
   storage_.addAtom();
   TEST_ASSERT(atomPtr22 == storage_.atomPtr(22));
   TEST_ASSERT(0 == storage_.atomPtr(15));
   TEST_ASSERT(0 == storage_.atomPtr(13));
   TEST_ASSERT(storage_.nAtom() == 1);

   Atom* atomPtr13;
   atomPtr13 = storage_.newAtomPtr();
   atomPtr13->id = 13;
   storage_.addAtom();
   TEST_ASSERT(atomPtr13 == storage_.atomPtr(13));
   TEST_ASSERT(atomPtr22 == storage_.atomPtr(22));
   TEST_ASSERT(0 == storage_.atomPtr(15));
   TEST_ASSERT(storage_.nAtom() == 2);

   Atom* atomPtr15;
   atomPtr15 = storage_.newAtomPtr();
   atomPtr15->id = 15;
   storage_.addAtom();
   TEST_ASSERT(atomPtr15 == storage_.atomPtr(15));
   TEST_ASSERT(atomPtr13 == storage_.atomPtr(13));
   TEST_ASSERT(atomPtr22 == storage_.atomPtr(22));
   TEST_ASSERT(storage_.nAtom() == 3);
}

TEST_BEGIN(StorageTest)
TEST_ADD(StorageTest, testReadParam)
TEST_ADD(StorageTest, testAddAtoms)
TEST_END(StorageTest)

#endif
