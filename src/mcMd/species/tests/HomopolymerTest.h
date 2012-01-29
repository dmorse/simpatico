#ifndef HOMOPOLYMER_TEST_H
#define HOMOPOLYMER_TEST_H

#include <test/UnitTest.h>

#include <mcMd/species/Homopolymer.h>

#include <iostream>
#include <fstream>

using namespace Util;
using namespace McMd;

class HomopolymerTest : public UnitTest 
{

private:

   Homopolymer  species;

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
   void testPopMolecule();

};



void HomopolymerTest::testConstructor()
{
   printMethod(TEST_FUNC);
} 

void HomopolymerTest::testReadParam()
{
   printMethod(TEST_FUNC);
   using std::ifstream;
   using std::cout;

   ifstream in("in/Homopolymer");
   species.readParam(in);
   in.close();

   if (verbose() > 1) {
      species.writeParam(std::cout);
   }

   TEST_ASSERT(species.isValid());

}


#if 0
void HomopolymerTest::testPopMolecule()
{
   printMethod(TEST_FUNC);
   using std::ifstream;
   using std::cout;

   //cout << "Species Id  = " << species.id() << std::endl;
   TEST_ASSERT(species.id() == speciesId);

   // Read input file  
   std::ifstream in("species/in/Homopolymer");
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

TEST_BEGIN(HomopolymerTest)
TEST_ADD(HomopolymerTest, testConstructor)
TEST_ADD(HomopolymerTest, testReadParam)
   //TEST_ADD(HomopolymerTest, testPopMolecule);
TEST_END(HomopolymerTest)

#endif
