#ifndef MDCF_STORAGE_TEST_H
#define MDCF_STORAGE_TEST_H

#include <mdCf/storage/System.h>
#include <mdCf/chemistry/Atom.h>
#include <mdCf/chemistry/Group.h>

#include <test/ParamFileTest.h>
#include <test/UnitTestRunner.h>

using namespace Util;
using namespace MdCf;

class SystemTest : public ParamFileTest
{

private:

   System system_;

public:

   SystemTest() 
    : system_()
   {}

   virtual void setUp()
   { 
      system_.readParam("in/System"); 
   }

   void testReadParam();

   void testAddAtoms();

};

inline void SystemTest::testReadParam()
{  
   printMethod(TEST_FUNC); 
   if (verbose() > 0) {
      system_.writeParam(std::cout);
   }
   TEST_ASSERT(system_.nSpecies() == 0);
}

inline void SystemTest::testAddAtoms()
{
   printMethod(TEST_FUNC);

   Atom* atomPtr22;
   atomPtr22 = system_.atoms().newPtr();
   atomPtr22->id = 22;
   system_.atoms().add();
   TEST_ASSERT(atomPtr22 == system_.atoms().ptr(22));
   TEST_ASSERT(0 == system_.atoms().ptr(15));
   TEST_ASSERT(0 == system_.atoms().ptr(13));
   TEST_ASSERT(system_.atoms().size() == 1);

   Atom* atomPtr13;
   atomPtr13 = system_.atoms().newPtr();
   atomPtr13->id = 13;
   system_.atoms().add();
   TEST_ASSERT(atomPtr13 == system_.atoms().ptr(13));
   TEST_ASSERT(atomPtr22 == system_.atoms().ptr(22));
   TEST_ASSERT(0 == system_.atoms().ptr(15));
   TEST_ASSERT(system_.atoms().size() == 2);

   Atom* atomPtr15;
   atomPtr15 = system_.atoms().newPtr();
   atomPtr15->id = 15;
   system_.atoms().add();
   TEST_ASSERT(atomPtr15 == system_.atoms().ptr(15));
   TEST_ASSERT(atomPtr13 == system_.atoms().ptr(13));
   TEST_ASSERT(atomPtr22 == system_.atoms().ptr(22));
   TEST_ASSERT(system_.atoms().size() == 3);
}

TEST_BEGIN(SystemTest)
TEST_ADD(SystemTest, testReadParam)
TEST_ADD(SystemTest, testAddAtoms)
TEST_END(SystemTest)

#endif
