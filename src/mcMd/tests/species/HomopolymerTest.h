#ifndef MCMD_HOMOPOLYMER_TEST_H
#define MCMD_HOMOPOLYMER_TEST_H

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

TEST_BEGIN(HomopolymerTest)
TEST_ADD(HomopolymerTest, testConstructor)
TEST_ADD(HomopolymerTest, testReadParam)
TEST_END(HomopolymerTest)

#endif
