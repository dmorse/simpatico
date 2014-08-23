#ifndef SPAN_PROCESSOR_TEST_H
#define SPAN_PROCESSOR_TEST_H

#include <spAn/processor/Processor.h>
#include <spAn/chemistry/Atom.h>
#include <spAn/chemistry/Group.h>

#include <test/ParamFileTest.h>
#include <test/UnitTestRunner.h>

using namespace Util;
using namespace SpAn;

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
   readParam("in/Processor");
   std::ifstream file;
   openInputFile("in/config", file);
   processor_.readConfig(file);
   file.close();
}

inline void ProcessorTest::testWriteConfig1()
{
   printMethod(TEST_FUNC);
   readParam("in/Processor");
   std::ofstream file;
   openOutputFile("out/config", file);
   processor_.writeConfig(file);
   file.close();
}

inline void ProcessorTest::testReadConfig2()
{
   printMethod(TEST_FUNC);
   //ParamComponent::setEcho(true);
   readParam("in/Processor.2");
   TEST_ASSERT(processor_.nSpecies() == 1);
   processor_.setConfigIo("DdMdConfigIo_Molecule");
   std::ifstream file;
   openInputFile("in/config.3", file);
   processor_.readConfig(file);
   file.close();
   TEST_ASSERT(processor_.species(0).size() == 8);
   TEST_ASSERT(processor_.atoms().size() == 256);

   #if 0
   Species* speciesPtr;
   Species::Iterator mIter; 
   Atom* atomPtr;
   int nSpecies = processor_.nSpecies();
   int i;
   for (int is = 0; is < nSpecies; ++is) {
      speciesPtr = &processor_.species(is);
      for (speciesPtr->begin(mIter); mIter.notEnd(); ++mIter) {
         for (i = 0; i < speciesPtr->nAtom(); ++i) {
            atomPtr = &mIter->atom(i);
            std::cout << atomPtr->id << " "
                      << atomPtr->typeId << " "
                      << atomPtr->speciesId << " "
                      << atomPtr->moleculeId << " "
                      << atomPtr->atomId << " "
                      << atomPtr->position << " "
                      << std::endl;
         }
      }
   }
   #endif
}

TEST_BEGIN(ProcessorTest)
TEST_ADD(ProcessorTest, testReadParam)
TEST_ADD(ProcessorTest, testAddAtoms)
TEST_ADD(ProcessorTest, testReadConfig1)
TEST_ADD(ProcessorTest, testWriteConfig1)
TEST_ADD(ProcessorTest, testReadConfig2)
TEST_END(ProcessorTest)

#endif
