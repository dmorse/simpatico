#ifndef MCMD_SPECIES_TEST_H
#define MCMD_SPECIES_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <simp/species/Species.h>

#include <iostream>
#include <fstream>

using namespace Util;
using namespace McMd;

class SpeciesTest : public UnitTest 
{

private:

   Species   species;

   const static int speciesId = 11; 

public:

   void setUp() 
   { 
      setVerbose(2);
      species.setId(speciesId); 
   } 

   void tearDown() 
   {}
  
   void testConstructor();

   #ifdef SIMP_BOND
   void testReadParamBond();
   #endif

   #ifdef SIMP_ANGLE
   void testReadParamAngle();
   #endif

   #ifdef SIMP_DIHEDRAL
   void testReadParamDihedral();
   #endif

   void testWriteStructure();

   void testReadStructure();

   void testMatchStructure();

};



void SpeciesTest::testConstructor()
{
   printMethod(TEST_FUNC);
} 

#ifdef SIMP_BOND
void SpeciesTest::testReadParamBond()
{
   printMethod(TEST_FUNC);

   std::ifstream in;
   openInputFile("in/Species", in);
   species.readParam(in);
   in.close();

   if (verbose() > 1) {
      std::cout << std::endl;
      species.writeParam(std::cout);
   }

   TEST_ASSERT(species.isValid());
}
#endif

#ifdef SIMP_ANGLE
void SpeciesTest::testReadParamAngle()
{
   printMethod(TEST_FUNC);

   std::ifstream in;
   openInputFile("in/SpeciesAngle", in);
   species.readParam(in);
   in.close();

   if (verbose() > 1) {
      std::cout << std::endl;
      species.writeParam(std::cout);
   }

   TEST_ASSERT(species.isValid());

}
#endif

#ifdef SIMP_DIHEDRAL
void SpeciesTest::testReadParamDihedral()
{
   printMethod(TEST_FUNC);

   std::ifstream in;
   openInputFile("in/SpeciesAngleDihedral", in);
   species.readParam(in);
   in.close();

   if (verbose() > 1) {
      std::cout << std::endl;
      species.writeParam(std::cout);
   }

   TEST_ASSERT(species.isValid());

}
#endif

void SpeciesTest::testWriteStructure()
{
   printMethod(TEST_FUNC);

   std::ifstream in;
   #ifdef SIMP_ANGLE
   #ifdef SIMP_DIHEDRAL
   openInputFile("in/SpeciesAngleDihedral", in);
   #else
   openInputFile("in/SpeciesAngle", in);
   #endif
   #else
   openInputFile("in/Species", in);
   #endif
   species.readParam(in);
   in.close();

   if (verbose() > 1) {
      std::cout << std::endl;
      species.writeParam(std::cout);
   }
   TEST_ASSERT(species.isValid());

   species.writeStructure(std::cout, "  ");

   std::ofstream out;
   openOutputFile("tmp/structure", out);
   species.writeStructure(out, "  ");

}

void SpeciesTest::testReadStructure()
{
   printMethod(TEST_FUNC);
   std::ifstream in;
   openInputFile("in/structure", in);
   species.readStructure(in);
   in.close();

   TEST_ASSERT(species.isValid());
   if (verbose() > 1) {
      std::cout << std::endl;
      species.writeStructure(std::cout, "  ");
   }
}

void SpeciesTest::testMatchStructure()
{
   printMethod(TEST_FUNC);
   std::ifstream in;
   openInputFile("in/structure", in);
   species.readStructure(in);
   in.close();

   TEST_ASSERT(species.isValid());

   openInputFile("in/structure", in);
   species.matchStructure(in);
   in.close();

}

TEST_BEGIN(SpeciesTest)
TEST_ADD(SpeciesTest, testConstructor)
#ifdef SIMP_BOND
TEST_ADD(SpeciesTest, testReadParamBond)
#endif
#ifdef SIMP_ANGLE
TEST_ADD(SpeciesTest, testReadParamAngle)
#endif
#ifdef SIMP_DIHEDRAL
TEST_ADD(SpeciesTest, testReadParamDihedral)
#endif
TEST_ADD(SpeciesTest, testWriteStructure)
TEST_ADD(SpeciesTest, testReadStructure)
TEST_ADD(SpeciesTest, testMatchStructure)
TEST_END(SpeciesTest)

#endif
