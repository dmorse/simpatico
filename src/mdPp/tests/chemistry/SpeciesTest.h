#ifndef MDPP_SPECIES_TEST_H
#define MDPP_SPECIES_TEST_H

#include <mdPp/chemistry/Species.h>
#include <mdPp/chemistry/Atom.h>
#include <mdPp/chemistry/Molecule.h>

#include <test/ParamFileTest.h>
#include <test/UnitTestRunner.h>

using namespace Util;
using namespace MdPp;

class SpeciesTest : public ParamFileTest
{

private:

   Species species_;

public:

   SpeciesTest() 
    : species_()
   {}

   virtual void setUp()
   {
      species_.initialize(3, 3); 
   }

   void testReadParam();

   void testAddAtoms();

};

inline void SpeciesTest::testAddAtoms()
{
   printMethod(TEST_FUNC);

   DArray<Atom> atoms;
   TEST_ASSERT(species_.nAtom() == 0);

}

TEST_BEGIN(SpeciesTest)
TEST_ADD(SpeciesTest, testAddAtoms)
TEST_END(SpeciesTest)

#endif
