#ifndef SPAN_HOOMD_CONFIG_READER_TEST_H
#define SPAN_HOOMD_CONFIG_READER_TEST_H

#include <spAn/processor/Processor.h>
#include <spAn/chemistry/Atom.h>
#include <spAn/chemistry/Group.h>

#include <test/ParamFileTest.h>
#include <test/UnitTestRunner.h>

using namespace Util;
using namespace SpAn;

class HoomdConfigReaderTest : public ParamFileTest
{

private:

   Processor processor_;

public:

   HoomdConfigReaderTest() 
    : processor_()
   {}

   virtual void setUp()
   { 
      setVerbose(2);
   }

   void testReadParam();
   void testReadConfig();

};

inline void HoomdConfigReaderTest::testReadParam()
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

inline void HoomdConfigReaderTest::testReadConfig()
{
   printMethod(TEST_FUNC);

   std::ifstream file;
   openInputFile("in/Processor", file);
   processor_.readParam(file); 
   file.close(); 
   
   processor_.setConfigReader("HoomdConfigReader");

   openInputFile("in/config.hoomd", file);
   processor_.readConfig(file);
   file.close(); 

   //TEST_ASSERT(processor_.atoms().size() == 40);
   //TEST_ASSERT(processor_.bonds().size() == 35);
   // std::cout << "nAtom = " << processor_.atoms().size() << "\n";
   // std::cout << "nBond = " << processor_.bonds().size() << "\n";
   // std::cout << processor_.boundary() << "\n";

}

TEST_BEGIN(HoomdConfigReaderTest)
TEST_ADD(HoomdConfigReaderTest, testReadParam)
TEST_ADD(HoomdConfigReaderTest, testReadConfig)
TEST_END(HoomdConfigReaderTest)

#endif
