#ifndef SPAN_HOOMD_CONFIG_READER_TEST_H
#define SPAN_HOOMD_CONFIG_READER_TEST_H

#include <spAn/processor/Processor.h>
#include <spAn/chemistry/Atom.h>
#include <spAn/chemistry/Group.h>

#include <test/ParamFileTest.h>
#include <test/UnitTestRunner.h>

#include <iostream>

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

   // Read auxiliary file
   openInputFile("in/types", file);
   processor_.configReader().readAuxiliaryFile(file);
   //std::cout << std::endl << "Read auxiliary types file";
   file.close(); 

   // Read auxiliary file
   openInputFile("in/config.hoomd", file);
   processor_.readConfig(file);
   file.close(); 

   TEST_ASSERT(processor_.atoms().size() == 32);
   TEST_ASSERT(processor_.bonds().size() == 31);
 
   #if 0
   std::cout << std::endl;
   std::cout << processor_.boundary() << "\n";

   AtomStorage::Iterator atomIt;
   std::cout << std::endl;
   for (processor_.atoms().begin(atomIt); atomIt.notEnd(); ++atomIt) {
      std::cout << atomIt->id << "  " 
                << atomIt->typeId << std::endl;
   }

   GroupStorage< 2 > ::Iterator bondIt;
   std::cout << std::endl;
   for (processor_.bonds().begin(bondIt); bondIt.notEnd(); ++bondIt) {
      std::cout << bondIt->id << "  " 
                << bondIt->typeId  << "  " 
                << bondIt->atomIds[0] << "  "
                << bondIt->atomIds[1] << "  "
                << std::endl;
   }
   #endif

   //  Open and write output hoomd file
   processor_.setConfigWriter("HoomdConfigWriter");
   openInputFile("in/types", file);
   processor_.configWriter().readAuxiliaryFile(file);
   file.close(); 
   std::ofstream out;
   openOutputFile("out/config.hoomd", out);
   processor_.writeConfig(out);
   out.close(); 

}

TEST_BEGIN(HoomdConfigReaderTest)
TEST_ADD(HoomdConfigReaderTest, testReadParam)
TEST_ADD(HoomdConfigReaderTest, testReadConfig)
TEST_END(HoomdConfigReaderTest)

#endif
