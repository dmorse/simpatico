#ifndef MDPP_PROCESSOR_TEST_H
#define MDPP_PROCESSOR_TEST_H

#include <mdPp/Processor.h>
//#include <mdPp/chemistry/Atom.h>
//#include <mdPp/chemistry/Group.h>
//#include <util/containers/DArray.h>

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

   virtual void setUp()
   { 
      processor_.readParam("in/Processor"); 
   }

   void testReadParam();

   void testAddAtoms();

};

inline void ProcessorTest::testReadParam()
{  
   printMethod(TEST_FUNC); 
   if (verbose() > 0) {
      processor_.writeParam(std::cout);
   }
}

inline void ProcessorTest::testAddAtoms()
{
   printMethod(TEST_FUNC);

   Atom* atomPtr22;
   atomPtr22 = processor_.newAtomPtr();
   atomPtr22->id = 22;
   processor_.addAtom();
   TEST_ASSERT(atomPtr22 == processor_.atomPtr(22));
   TEST_ASSERT(0 == processor_.atomPtr(15));
   TEST_ASSERT(0 == processor_.atomPtr(13));
   TEST_ASSERT(processor_.nAtom() == 1);
   TEST_ASSERT(processor_.nBond() == 0);

   Atom* atomPtr13;
   atomPtr13 = processor_.newAtomPtr();
   atomPtr13->id = 13;
   processor_.addAtom();
   TEST_ASSERT(atomPtr13 == processor_.atomPtr(13));
   TEST_ASSERT(atomPtr22 == processor_.atomPtr(22));
   TEST_ASSERT(0 == processor_.atomPtr(15));
   TEST_ASSERT(processor_.nAtom() == 2);

   Atom* atomPtr15;
   atomPtr15 = processor_.newAtomPtr();
   atomPtr15->id = 15;
   processor_.addAtom();
   TEST_ASSERT(atomPtr15 == processor_.atomPtr(15));
   TEST_ASSERT(atomPtr13 == processor_.atomPtr(13));
   TEST_ASSERT(atomPtr22 == processor_.atomPtr(22));
   TEST_ASSERT(processor_.nAtom() == 3);
   TEST_ASSERT(processor_.nBond() == 0);
}

TEST_BEGIN(ProcessorTest)
TEST_ADD(ProcessorTest, testReadParam)
TEST_ADD(ProcessorTest, testAddAtoms)
TEST_END(ProcessorTest)

#endif
