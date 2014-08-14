#ifndef SPAN_DDMD_CONFIGIO_TEST_H
#define SPAN_DDMD_CONFIGIO_TEST_H

#include <spAn/processor/Processor.h>
#include <spAn/chemistry/Atom.h>
#include <spAn/chemistry/Group.h>

#include <test/ParamFileTest.h>
#include <test/UnitTestRunner.h>

using namespace Util;
using namespace SpAn;

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
   std::ifstream file;
   openInputFile("in/Processor", file);
   processor_.readParam(file); 
   file.close(); 
   if (verbose() > 0) {
      std::cout << "\n";
      processor_.writeParam(std::cout);
   }
}

inline void DdMdConfigIoTest::testReadConfig1()
{
   printMethod(TEST_FUNC);

   std::ifstream file;
   openInputFile("in/Processor", file);
   processor_.readParam(file); 
   file.close(); 

   openInputFile("in/config", file);
   processor_.readConfig(file);
   file.close(); 

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
   std::ifstream file;

   openInputFile("in/Processor.2", file);
   processor_.readParam(file); 
   file.close();

   processor_.setConfigIo("DdMdConfigIo_Molecule");
   openInputFile("in/config.2", file);
   processor_.readConfig(file);
   file.close();

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
