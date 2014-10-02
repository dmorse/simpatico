#ifndef SPAN_STORAGE_TEST_H
#define SPAN_STORAGE_TEST_H

#include <spAn/storage/Configuration.h>
#include <spAn/chemistry/Atom.h>
#include <spAn/chemistry/Group.h>

#include <test/ParamFileTest.h>
#include <test/UnitTestRunner.h>

using namespace Util;
using namespace SpAn; 

class ConfigurationTest : public ParamFileTest
{

private:

   Configuration configuration_;

public:

   ConfigurationTest() 
    : configuration_()
   {}

   virtual void setUp()
   { 
      std::ifstream file;
      openInputFile("in/Configuration", file);
      configuration_.readParam(file); 
      file.close(); 
      //configuration_.readParam("in/Configuration"); 
   }

   void testReadParam();

   void testAddAtoms();

};

inline void ConfigurationTest::testReadParam()
{  
   printMethod(TEST_FUNC); 
   if (verbose() > 0) {
      configuration_.writeParam(std::cout);
   }
   TEST_ASSERT(configuration_.nSpecies() == 0);
}

inline void ConfigurationTest::testAddAtoms()
{
   printMethod(TEST_FUNC);

   Atom* atomPtr22;
   atomPtr22 = configuration_.atoms().newPtr();
   atomPtr22->id = 22;
   configuration_.atoms().add();
   TEST_ASSERT(atomPtr22 == configuration_.atoms().ptr(22));
   TEST_ASSERT(0 == configuration_.atoms().ptr(15));
   TEST_ASSERT(0 == configuration_.atoms().ptr(13));
   TEST_ASSERT(configuration_.atoms().size() == 1);

   Atom* atomPtr13;
   atomPtr13 = configuration_.atoms().newPtr();
   atomPtr13->id = 13;
   configuration_.atoms().add();
   TEST_ASSERT(atomPtr13 == configuration_.atoms().ptr(13));
   TEST_ASSERT(atomPtr22 == configuration_.atoms().ptr(22));
   TEST_ASSERT(0 == configuration_.atoms().ptr(15));
   TEST_ASSERT(configuration_.atoms().size() == 2);

   Atom* atomPtr15;
   atomPtr15 = configuration_.atoms().newPtr();
   atomPtr15->id = 15;
   configuration_.atoms().add();
   TEST_ASSERT(atomPtr15 == configuration_.atoms().ptr(15));
   TEST_ASSERT(atomPtr13 == configuration_.atoms().ptr(13));
   TEST_ASSERT(atomPtr22 == configuration_.atoms().ptr(22));
   TEST_ASSERT(configuration_.atoms().size() == 3);
}

TEST_BEGIN(ConfigurationTest)
TEST_ADD(ConfigurationTest, testReadParam)
TEST_ADD(ConfigurationTest, testAddAtoms)
TEST_END(ConfigurationTest)

#endif
