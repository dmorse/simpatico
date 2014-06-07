#ifndef DDMD_SP_STORAGE_TEST_H
#define DDMD_SP_STORAGE_TEST_H

#include <ddMd/sp/storage/SpConfiguration.h>
#include <ddMd/sp/chemistry/SpAtom.h>
#include <ddMd/sp/chemistry/SpGroup.h>

#include <test/ParamFileTest.h>
#include <test/UnitTestRunner.h>

using namespace Util;
using namespace DdMd; 

class SpConfigurationTest : public ParamFileTest
{

private:

   SpConfiguration configuration_;

public:

   SpConfigurationTest() 
    : configuration_()
   {}

   virtual void setUp()
   { 
      std::ifstream file;
      openInputFile("in/SpConfiguration", file);
      configuration_.readParam(file); 
      file.close(); 
      //configuration_.readParam("in/SpConfiguration"); 
   }

   void testReadParam();

   void testAddAtoms();

};

inline void SpConfigurationTest::testReadParam()
{  
   printMethod(TEST_FUNC); 
   if (verbose() > 0) {
      configuration_.writeParam(std::cout);
   }
   TEST_ASSERT(configuration_.nSpecies() == 0);
}

inline void SpConfigurationTest::testAddAtoms()
{
   printMethod(TEST_FUNC);

   SpAtom* atomPtr22;
   atomPtr22 = configuration_.atoms().newPtr();
   atomPtr22->id = 22;
   configuration_.atoms().add();
   TEST_ASSERT(atomPtr22 == configuration_.atoms().ptr(22));
   TEST_ASSERT(0 == configuration_.atoms().ptr(15));
   TEST_ASSERT(0 == configuration_.atoms().ptr(13));
   TEST_ASSERT(configuration_.atoms().size() == 1);

   SpAtom* atomPtr13;
   atomPtr13 = configuration_.atoms().newPtr();
   atomPtr13->id = 13;
   configuration_.atoms().add();
   TEST_ASSERT(atomPtr13 == configuration_.atoms().ptr(13));
   TEST_ASSERT(atomPtr22 == configuration_.atoms().ptr(22));
   TEST_ASSERT(0 == configuration_.atoms().ptr(15));
   TEST_ASSERT(configuration_.atoms().size() == 2);

   SpAtom* atomPtr15;
   atomPtr15 = configuration_.atoms().newPtr();
   atomPtr15->id = 15;
   configuration_.atoms().add();
   TEST_ASSERT(atomPtr15 == configuration_.atoms().ptr(15));
   TEST_ASSERT(atomPtr13 == configuration_.atoms().ptr(13));
   TEST_ASSERT(atomPtr22 == configuration_.atoms().ptr(22));
   TEST_ASSERT(configuration_.atoms().size() == 3);
}

TEST_BEGIN(SpConfigurationTest)
TEST_ADD(SpConfigurationTest, testReadParam)
TEST_ADD(SpConfigurationTest, testAddAtoms)
TEST_END(SpConfigurationTest)

#endif
