#ifndef MCMD_SPECIES_GROUP_TEST_H
#define MCMD_SPECIES_GROUP_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <simp/species/SpeciesGroup.tpp>
#include <util/archives/BinaryFileOArchive.h>
#include <util/archives/BinaryFileIArchive.h>

using namespace Simp;

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

   void testSerialize() {
      printMethod(TEST_FUNC);
      SpeciesGroup<2> v;
      int i1 = 35;
      int i2 = 43;

      // Read from input file
      std::ifstream in;
      openInputFile("in/SpeciesGroup", in);
      in        >> v;

      // Write to binary file archive
      BinaryFileOArchive oa;
      openOutputFile("binary", oa.file());
      oa << i1;
      oa << v;
      oa << i2;
      oa.file().close();

      // Write to binary file archive
      SpeciesGroup<2> u;
      int j1, j2;
      BinaryFileIArchive ia;
      openInputFile("binary", ia.file());
      ia >> j1;
      ia >> u;
      ia >> j2;
      
      TEST_ASSERT(j1 == i1);
      TEST_ASSERT(u.typeId() == v.typeId());
      //std::cout << v.typeId() << std::endl;
      for (int i = 0; i < 2; ++i) {
         TEST_ASSERT(u.atomId(i) == v.atomId(i));
         //std::cout << v.atomId(i) << std::endl;
      }
      TEST_ASSERT(j2 == i2);
   }

};

TEST_BEGIN(SpeciesGroupTest)
TEST_ADD(SpeciesGroupTest, testReadWrite)
TEST_ADD(SpeciesGroupTest, testSerialize)
TEST_END(SpeciesGroupTest)

#endif
