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
   //void testPopMolecule();

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

#if 0
void SpeciesTest::testPopMolecule()
{
   printMethod(TEST_FUNC);
   using std::ifstream;
   using std::cout;

   //cout << "Species Id  = " << species.id() << std::endl;
   TEST_ASSERT(species.id() == speciesId);

   // Read input file  
   std::ifstream in("species/in/Species");
   species.readParam(in);
   in.close();

   TEST_ASSERT(species.isValid());

   Molecule* molPtr;

   molPtr = &(species.reservoir().pop());
   //cout << "Molecule Id = " << molPtr->id() << std::endl;
   //cout << "Species Id  = " << molPtr->species().id() << std::endl;
   TEST_ASSERT(&(molPtr->species()) == &species );
   TEST_ASSERT(molPtr->id() == 0 );
   TEST_ASSERT(species.reservoir().size() 
                   == species.reservoir().capacity() - 1);
   TEST_ASSERT(species.isValid() );

   molPtr = &(species.reservoir().pop());
   //cout << "Molecule Id = " << molPtr->id() << std::endl;
   //cout << "Species Id  = " << molPtr->species().id() << std::endl;
   TEST_ASSERT(&(molPtr->species()) == &species);
   TEST_ASSERT(molPtr->id() == 1);
   TEST_ASSERT(species.reservoir().size() 
                   == species.reservoir().capacity() - 2);

   molPtr = &(species.reservoir().pop());
   //cout << "Molecule Id = " << molPtr->id() << std::endl;
   //cout << "Species Id  = " << molPtr->species().id() << std::endl;
   TEST_ASSERT(&(molPtr->species()) == &species);
   TEST_ASSERT(molPtr->id() == 2);
   TEST_ASSERT(species.reservoir().size() 
                   == species.reservoir().capacity() - 3);

   molPtr = &(species.reservoir().pop());
   //cout << "Molecule Id = " << molPtr->id() << std::endl;
   //cout << "Species Id  = " << molPtr->species().id() << std::endl;
   TEST_ASSERT(&(molPtr->species()) == &species);
   TEST_ASSERT(molPtr->id() == 3);
   TEST_ASSERT(species.reservoir().size() 
                   == species.reservoir().capacity() - 4);

   TEST_ASSERT(species.isValid());

}
#endif

TEST_BEGIN(SpeciesTest)
TEST_ADD(SpeciesTest, testConstructor)
TEST_ADD(SpeciesTest, testReadParam)
   //TEST_ADD(SpeciesTest, testPopMolecule);
TEST_END(SpeciesTest)

#endif
