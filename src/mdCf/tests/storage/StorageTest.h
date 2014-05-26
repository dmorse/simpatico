#ifndef MDPP_STORAGE_TEST_H
#define MDPP_STORAGE_TEST_H

#include <mdCf/storage/Storage.h>
#include <mdCf/chemistry/Atom.h>
#include <mdCf/chemistry/Group.h>

#include <test/ParamFileTest.h>
#include <test/UnitTestRunner.h>

using namespace Util;
using namespace MdCf;

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
   atomPtr22 = storage_.atoms().newPtr();
   atomPtr22->id = 22;
   storage_.atoms().add();
   TEST_ASSERT(atomPtr22 == storage_.atoms().ptr(22));
   TEST_ASSERT(0 == storage_.atoms().ptr(15));
   TEST_ASSERT(0 == storage_.atoms().ptr(13));
   TEST_ASSERT(storage_.atoms().size() == 1);

   Atom* atomPtr13;
   atomPtr13 = storage_.atoms().newPtr();
   atomPtr13->id = 13;
   storage_.atoms().add();
   TEST_ASSERT(atomPtr13 == storage_.atoms().ptr(13));
   TEST_ASSERT(atomPtr22 == storage_.atoms().ptr(22));
   TEST_ASSERT(0 == storage_.atoms().ptr(15));
   TEST_ASSERT(storage_.atoms().size() == 2);

   Atom* atomPtr15;
   atomPtr15 = storage_.atoms().newPtr();
   atomPtr15->id = 15;
   storage_.atoms().add();
   TEST_ASSERT(atomPtr15 == storage_.atoms().ptr(15));
   TEST_ASSERT(atomPtr13 == storage_.atoms().ptr(13));
   TEST_ASSERT(atomPtr22 == storage_.atoms().ptr(22));
   TEST_ASSERT(storage_.atoms().size() == 3);
}

TEST_BEGIN(StorageTest)
TEST_ADD(StorageTest, testReadParam)
TEST_ADD(StorageTest, testAddAtoms)
TEST_END(StorageTest)

#endif
