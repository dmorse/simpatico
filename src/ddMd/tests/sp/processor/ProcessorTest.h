#ifndef DDMD_SP_PROCESSOR_TEST_H
#define DDMD_SP_PROCESSOR_TEST_H

#include <ddMd/sp/processor/Processor.h>
#include <ddMd/sp/chemistry/SpAtom.h>
#include <ddMd/sp/chemistry/SpGroup.h>

#include <test/ParamFileTest.h>
#include <test/UnitTestRunner.h>

using namespace Util;
using namespace DdMd;

class ProcessorTest : public ParamFileTest
{

private:

   Processor processor_;

public:

   ProcessorTest() 
    : processor_()
   {}

   virtual void setUp()
   { 
      std::ifstream file;
      openInputFile("in/Processor", file);
      processor_.readParam(file); 
      file.close(); 
   }

   void testReadParam();
   void testAddAtoms();
   void testReadConfig();

};

inline void ProcessorTest::testReadParam()
{  
   printMethod(TEST_FUNC); 
   if (verbose() > 0) {
      processor_.writeParam(std::cout);
   }
   TEST_ASSERT(processor_.nSpecies() == 0);
}

inline void ProcessorTest::testAddAtoms()
{
   printMethod(TEST_FUNC);

   SpAtom* atomPtr22;
   atomPtr22 = processor_.atoms().newPtr();
   atomPtr22->id = 22;
   processor_.atoms().add();
   TEST_ASSERT(atomPtr22 == processor_.atoms().ptr(22));
   TEST_ASSERT(0 == processor_.atoms().ptr(15));
   TEST_ASSERT(0 == processor_.atoms().ptr(13));
   TEST_ASSERT(processor_.atoms().size() == 1);

   SpAtom* atomPtr13;
   atomPtr13 = processor_.atoms().newPtr();
   atomPtr13->id = 13;
   processor_.atoms().add();
   TEST_ASSERT(atomPtr13 == processor_.atoms().ptr(13));
   TEST_ASSERT(atomPtr22 == processor_.atoms().ptr(22));
   TEST_ASSERT(0 == processor_.atoms().ptr(15));
   TEST_ASSERT(processor_.atoms().size() == 2);

   SpAtom* atomPtr15;
   atomPtr15 = processor_.atoms().newPtr();
   atomPtr15->id = 15;
   processor_.atoms().add();
   TEST_ASSERT(atomPtr15 == processor_.atoms().ptr(15));
   TEST_ASSERT(atomPtr13 == processor_.atoms().ptr(13));
   TEST_ASSERT(atomPtr22 == processor_.atoms().ptr(22));
   TEST_ASSERT(processor_.atoms().size() == 3);
}

inline void ProcessorTest::testReadConfig()
{
  printMethod(TEST_FUNC);
  processor_.readConfig("in/config");
  processor_.writeConfig("out/config.out");
}

TEST_BEGIN(ProcessorTest)
TEST_ADD(ProcessorTest, testReadParam)
TEST_ADD(ProcessorTest, testAddAtoms)
TEST_ADD(ProcessorTest, testReadConfig)
TEST_END(ProcessorTest)

#endif
