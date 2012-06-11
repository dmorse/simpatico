#ifndef MCMD_ATOM_TYPE_TEST_H
#define MCMD_ATOM_TYPE_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <mcMd/chemistry/AtomType.h>

#include <fstream>

using namespace McMd;

class AtomTypeTest : public UnitTest 
{

public:

   void setUp()
   {}

   void tearDown()
   {}

   void testReadWrite() {
      printMethod(TEST_FUNC);
      AtomType v;
      std::ifstream in;
      openInputFile("in/AtomType", in);
      in        >> v;

      std::cout << std::endl;
      std::cout << v << std::endl ;
   }

};

TEST_BEGIN(AtomTypeTest )
TEST_ADD(AtomTypeTest, testReadWrite )
TEST_END(AtomTypeTest)


#endif
