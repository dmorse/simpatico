#ifndef MCMD_MASK_POLICY_TEST_H
#define MCMD_MASK_POLICY_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <mcMd/chemistry/MaskPolicy.h>
#include <util/archives/TextFileOArchive.h>
#include <util/archives/TextFileIArchive.h>

#include <fstream>

using namespace McMd;

class MaskPolicyTest : public UnitTest 
{


public:

   void setUp()
   {}

   void tearDown()
   {}

   void testReadWrite() {
      printMethod(TEST_FUNC);
      MaskPolicy policy;
      std::ifstream in;
      openInputFile("in/MaskPolicy", in);
      in >> policy;

      std::cout << std::endl;
      std::cout << policy << std::endl ;

      Util::TextFileOArchive oar;
      openOutputFile("ar.txt", oar.file());
      oar & policy;
      oar.file().close();

      MaskPolicy clone;
      Util::TextFileIArchive iar;
      openInputFile("ar.txt", iar.file());
      iar & clone;

      TEST_ASSERT(policy == clone);
      std::cout << clone;
   }

};

TEST_BEGIN(MaskPolicyTest)
TEST_ADD(MaskPolicyTest, testReadWrite )
TEST_END(MaskPolicyTest)

#endif
