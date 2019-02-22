#ifndef MDPP_SMP_CONFIG_READER_TEST_H
#define MDPP_SMP_CONFIG_READER_TEST_H

#include <mdPp/processor/Processor.h>
#include <mdPp/chemistry/Atom.h>
#include <mdPp/chemistry/Group.h>

#include <test/ParamFileTest.h>
#include <test/UnitTestRunner.h>

using namespace Util;
using namespace MdPp;

class SmpConfigReaderTest : public ParamFileTest
{

private:

   Processor processor_;

public:

   SmpConfigReaderTest() 
    : processor_()
   {}

   virtual void setUp()
   { 
      setVerbose(2);
   }

   void testReadConfig1();
   void testReadConfig2();

};

inline void SmpConfigReaderTest::testReadConfig1()
{
   printMethod(TEST_FUNC);

   std::ifstream file;
   //openInputFile("in/Processor", file);
   //processor_.readParam(file); 
   //file.close(); 

   // Read Smp (simpatico) configuration file 
   processor_.setConfigReader("SmpConfigReader");
   openInputFile("in/config.smc", file);
   processor_.readConfig(file);
   file.close(); 

   TEST_ASSERT(processor_.atoms().size() == 800);
   TEST_ASSERT(processor_.bonds().size() == 700);
   // std::cout << "nAtom = " << processor_.atoms().size() << "\n";
   // std::cout << "nBond = " << processor_.bonds().size() << "\n";
   // std::cout << processor_.boundary() << "\n";

   // Write in DdMd format
   processor_.setConfigWriter("DdMdConfigWriter");
   std::ofstream out;
   openOutputFile("out/config.ddc", out);
   processor_.writeConfig(out);
   out.close(); 

}

TEST_BEGIN(SmpConfigReaderTest)
TEST_ADD(SmpConfigReaderTest, testReadConfig1)
TEST_END(SmpConfigReaderTest)

#endif
