#ifndef SPECIES_GROUP_TEST_H
#define SPECIES_GROUP_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <mcMd/chemistry/SpeciesGroup.tpp>

using namespace McMd;

class SpeciesGroupTest : public UnitTest 
{

public:

   void setUp()
   {}

   void tearDown()
   {}

   void testReadWrite() {
      printMethod(TEST_FUNC);
      SpeciesGroup<2> v;
      std::ifstream in;
      openInputFile("in/SpeciesGroup", in);
      in        >> v;

      std::cout << std::endl ;
      std::cout << v << std::endl ;
   }

};

TEST_BEGIN(SpeciesGroupTest)
TEST_ADD(SpeciesGroupTest, testReadWrite)
TEST_END(SpeciesGroupTest)

#endif
