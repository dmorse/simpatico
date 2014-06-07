#ifndef DDMD_SP_DDMD_CONFIGIO_TEST_H
#define DDMD_SP_DDMD_CONFIGIO_TEST_H

#include <ddMd/sp/processor/Processor.h>
#include <ddMd/sp/chemistry/SpAtom.h>
#include <ddMd/sp/chemistry/SpGroup.h>

#include <test/ParamFileTest.h>
#include <test/UnitTestRunner.h>

using namespace Util;
using namespace DdMd;

class DdMdSpConfigIoTest : public ParamFileTest
{

private:

   Processor processor_;

public:

   DdMdSpConfigIoTest() 
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

inline void DdMdSpConfigIoTest::testReadParam()
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

inline void DdMdSpConfigIoTest::testReadConfig1()
{
   printMethod(TEST_FUNC);

   std::ifstream file;
   openInputFile("in/Processor", file);
   processor_.readParam(file); 
   file.close(); 

   //processor_.setSpConfigIo(std::string("DdMdSpConfigIo"));
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

inline void DdMdSpConfigIoTest::testReadConfig2()
{
   printMethod(TEST_FUNC);

   processor_.readParam("in/Processor.2"); 
   //processor_.setSpConfigIo(std::string("DdMdSpConfigIo_Molecule"));
   processor_.setSpConfigIo("DdMdSpConfigIo_Molecule");
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


TEST_BEGIN(DdMdSpConfigIoTest)
TEST_ADD(DdMdSpConfigIoTest, testReadParam)
TEST_ADD(DdMdSpConfigIoTest, testReadConfig1)
TEST_ADD(DdMdSpConfigIoTest, testReadConfig2)
TEST_END(DdMdSpConfigIoTest)

#endif
