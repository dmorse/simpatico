#ifndef MDPP_PROCESSOR_TEST_H
#define MDPP_PROCESSOR_TEST_H

#include <mdPp/processor/Processor.h>
#include <mdPp/chemistry/Atom.h>
#include <mdPp/chemistry/Group.h>

#include <test/ParamFileTest.h>
#include <test/UnitTestRunner.h>

using namespace Util;
using namespace MdPp;

class ProcessorTest : public ParamFileTest
{

private:

   Processor processor_;

public:

   ProcessorTest() 
    : processor_()
   {}

   void readParam(const char* filename);

   void testReadParam();
   void testAddAtoms();
   void testReadConfig1();
   void testWriteConfig1();
   void testReadConfig2();

};

inline void ProcessorTest::readParam(const char* filename)
{ 
   std::ifstream file;
   openInputFile(filename, file);
   processor_.readParam(file); 
   file.close(); 
}

inline void ProcessorTest::testReadParam()
{
   printMethod(TEST_FUNC); 
   readParam("in/Processor");
   if (verbose() > 0) {
      processor_.writeParam(std::cout);
   }
   TEST_ASSERT(processor_.nSpecies() == 0);
}

inline void ProcessorTest::testAddAtoms()
{
   printMethod(TEST_FUNC);
   readParam("in/Processor");

   Atom* atomPtr22;
   atomPtr22 = processor_.atoms().newPtr();
   atomPtr22->id = 22;
   processor_.atoms().add();
   TEST_ASSERT(atomPtr22 == processor_.atoms().ptr(22));
   TEST_ASSERT(0 == processor_.atoms().ptr(15));
   TEST_ASSERT(0 == processor_.atoms().ptr(13));
   TEST_ASSERT(processor_.atoms().size() == 1);

   Atom* atomPtr13;
   atomPtr13 = processor_.atoms().newPtr();
   atomPtr13->id = 13;
   processor_.atoms().add();
   TEST_ASSERT(atomPtr13 == processor_.atoms().ptr(13));
   TEST_ASSERT(atomPtr22 == processor_.atoms().ptr(22));
   TEST_ASSERT(0 == processor_.atoms().ptr(15));
   TEST_ASSERT(processor_.atoms().size() == 2);

   Atom* atomPtr15;
   atomPtr15 = processor_.atoms().newPtr();
   atomPtr15->id = 15;
   processor_.atoms().add();
   TEST_ASSERT(atomPtr15 == processor_.atoms().ptr(15));
   TEST_ASSERT(atomPtr13 == processor_.atoms().ptr(13));
   TEST_ASSERT(atomPtr22 == processor_.atoms().ptr(22));
   TEST_ASSERT(processor_.atoms().size() == 3);
}

inline void ProcessorTest::testReadConfig1()
{
   printMethod(TEST_FUNC);
   //ParamComponent::setEcho(true);
   //readParam("in/Processor.2");
   std::ifstream file;
   openInputFile("in/config", file);
   processor_.readConfig(file);
   file.close();
}

inline void ProcessorTest::testReadConfig2()
{
   printMethod(TEST_FUNC);
   processor_.setConfigReader("DdMdConfigReader_Molecule");
   std::ifstream file;
   openInputFile("in/config.3", file);
   processor_.readConfig(file);
   file.close();
   TEST_ASSERT(processor_.atoms().size() == 256);

}

TEST_BEGIN(ProcessorTest)
TEST_ADD(ProcessorTest, testReadParam)
TEST_ADD(ProcessorTest, testAddAtoms)
TEST_ADD(ProcessorTest, testReadConfig1)
TEST_ADD(ProcessorTest, testReadConfig2)
TEST_END(ProcessorTest)

#endif
