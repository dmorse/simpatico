#ifndef MCMD_HOMORING_TEST_H
#define MCMD_HOMORING_TEST_H

#include <test/UnitTest.h>

#include <simp/species/HomoRing.h>

#include <iostream>
#include <fstream>

using namespace Util;
using namespace Simp;

class HomoRingTest : public UnitTest 
{

private:

   HomoRing  species;

   const static int speciesId = 11;

public:

   void setUp() 
   { 
      species.setId(speciesId); 
      // setVerbose(2);
   } 

   void tearDown() 
   {}
  
   void testConstructor();
   void testReadParam();
   void testPopMolecule();

};



void HomoRingTest::testConstructor()
{
   printMethod(TEST_FUNC);
} 

void HomoRingTest::testReadParam()
{
   printMethod(TEST_FUNC);

   //ifstream in("in/HomoRing");
   std::ifstream in;
   openInputFile("in/HomoRing", in);
   species.readParam(in);
   in.close();

   if (verbose() > 1) {
      std::cout << std::endl;
      species.writeParam(std::cout);
   }

   TEST_ASSERT(species.isValid());

}

TEST_BEGIN(HomoRingTest)
TEST_ADD(HomoRingTest, testConstructor)
TEST_ADD(HomoRingTest, testReadParam)
TEST_END(HomoRingTest)

#endif
