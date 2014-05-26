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
      setVerbose(2);
   }

   void testReadParam();
   void testReadConfig1();
   void testReadConfig2();

};

inline void DdMdConfigIoTest::testReadParam()
{  
   printMethod(TEST_FUNC); 
   processor_.readParam("in/Processor"); 
   if (verbose() > 0) {
      std::cout << "\n";
      processor_.writeParam(std::cout);
   }
}

inline void DdMdConfigIoTest::testReadConfig1()
{
   printMethod(TEST_FUNC);

   processor_.readParam("in/Processor"); 

   //processor_.setConfigIo(std::string("DdMdConfigIo"));
   processor_.readConfig("in/config");

   TEST_ASSERT(processor_.atoms().size() == 40);
   TEST_ASSERT(processor_.bonds().size() == 35);
   // std::cout << "nAtom = " << processor_.atoms().size() << "\n";
   // std::cout << "nBond = " << processor_.bonds().size() << "\n";
   // std::cout << processor_.boundary() << "\n";

   std::ofstream out;
   openOutputFile("out/config", out);
   processor_.writeConfig(out);
   out.close();
}

inline void DdMdConfigIoTest::testReadConfig2()
{
   printMethod(TEST_FUNC);

   processor_.readParam("in/Processor.2"); 
   //processor_.setConfigIo(std::string("DdMdConfigIo_Molecule"));
   processor_.setConfigIo("DdMdConfigIo_Molecule");
   processor_.readConfig("in/config.2");

   TEST_ASSERT(processor_.atoms().size() == 256);
   TEST_ASSERT(processor_.bonds().size() == 248);
   // std::cout << "nAtom = " << processor_.atoms().size() << "\n";
   // std::cout << "nBond = " << processor_.bonds().size() << "\n";
   // std::cout << processor_.boundary() << "\n";

   std::ofstream out;
   openOutputFile("out/config.2", out);
   processor_.writeConfig(out);
   out.close();
}


TEST_BEGIN(DdMdConfigIoTest)
TEST_ADD(DdMdConfigIoTest, testReadParam)
TEST_ADD(DdMdConfigIoTest, testReadConfig1)
TEST_ADD(DdMdConfigIoTest, testReadConfig2)
TEST_END(DdMdConfigIoTest)

#endif
