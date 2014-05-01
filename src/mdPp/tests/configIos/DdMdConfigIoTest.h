#ifndef MDPP_DDMD_CONFIGIO_TEST_H
#define MDPP_DDMD_CONFIGIO_TEST_H

#include <mdPp/processor/Processor.h>
#include <mdPp/chemistry/Atom.h>
#include <mdPp/chemistry/Group.h>

#include <test/ParamFileTest.h>
#include <test/UnitTestRunner.h>

using namespace Util;
using namespace MdPp;

class DdMdConfigIoTest : public ParamFileTest
{

private:

   Processor processor_;

public:

   DdMdConfigIoTest() 
    : processor_()
   {}

   virtual void setUp()
   { 
      processor_.readParam("in/Processor"); 
      setVerbose(2);
   }

   void testReadParam();
   void testReadConfig();

};

inline void DdMdConfigIoTest::testReadParam()
{  
   printMethod(TEST_FUNC); 
   if (verbose() > 0) {
      std::cout << "\n";
      processor_.writeParam(std::cout);
   }
}

inline void DdMdConfigIoTest::testReadConfig()
{
   printMethod(TEST_FUNC);

   //processor_.setConfigIo(std::string("DdMdConfigIo"));

   processor_.readConfig("in/config");

   TEST_ASSERT(processor_.nAtom() == 40);
   TEST_ASSERT(processor_.bonds().size() == 35);
   // std::cout << "nAtom = " << processor_.nAtom() << "\n";
   // std::cout << "nBond = " << processor_.bonds().size() << "\n";
   // std::cout << processor_.boundary() << "\n";

   std::ofstream out;
   openOutputFile("out/config", out);
   processor_.writeConfig(out);
   out.close();
}

TEST_BEGIN(DdMdConfigIoTest)
TEST_ADD(DdMdConfigIoTest, testReadParam)
TEST_ADD(DdMdConfigIoTest, testReadConfig)
TEST_END(DdMdConfigIoTest)

#endif
