#ifndef MCMD_SPECIES_TEST_H
#define MCMD_SPECIES_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <mcMd/species/Species.h>

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
      species.setId(speciesId); 
   } 

   void tearDown() 
   {}
  
   void testConstructor();
   void testReadParam();

};



void SpeciesTest::testConstructor()
{
   printMethod(TEST_FUNC);
} 

void SpeciesTest::testReadParam()
{
   printMethod(TEST_FUNC);
   using std::ifstream;
   using std::cout;

   #ifdef INTER_ANGLE
   #ifdef INTER_DIHEDRAL
   ifstream in("in/SpeciesAngleDihedral");
   #else
   ifstream in("in/SpeciesAngle");
   #endif
   #else
   ifstream in("in/Species");
   #endif
   species.readParam(in);
   in.close();

   if (verbose() > 1) {
      species.writeParam(std::cout);
   }

   TEST_ASSERT(species.isValid());

}

TEST_BEGIN(SpeciesTest)
TEST_ADD(SpeciesTest, testConstructor)
TEST_ADD(SpeciesTest, testReadParam)
TEST_END(SpeciesTest)

#endif
