#ifndef MCMD_MULTIBLOCK_TEST_H
#define MCMD_MULTIBLOCK_TEST_H

#include <test/UnitTest.h>

#include <simp/species/Multiblock.h>

#include <iostream>
#include <fstream>

using namespace Util;
using namespace Simp;

class MultiblockTest : public UnitTest 
{

private:

   Multiblock  species;

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

};

void MultiblockTest::testConstructor()
{
   printMethod(TEST_FUNC);
} 

void MultiblockTest::testReadParam()
{
   printMethod(TEST_FUNC);
   using std::ifstream;
   using std::cout;

   ifstream in;
   #ifndef SIMP_ANGLE
   openInputFile("in/Multiblock", in);
   #else
   openInputFile("in/MultiblockAngle", in);
   #endif
   species.readParam(in);
   in.close();

   if (verbose() > 1) {
      std::cout << std::endl;
      species.writeParam(std::cout);
   }

   TEST_ASSERT(species.isValid());

   // Atom types
   TEST_ASSERT(species.nAtom() == 12);
   TEST_ASSERT(species.atomTypeId(0) == 0);
   TEST_ASSERT(species.atomTypeId(1) == 0);
   TEST_ASSERT(species.atomTypeId(2) == 0);
   TEST_ASSERT(species.atomTypeId(3) == 0);
   TEST_ASSERT(species.atomTypeId(4) == 1);
   TEST_ASSERT(species.atomTypeId(5) == 1);
   TEST_ASSERT(species.atomTypeId(6) == 1);
   TEST_ASSERT(species.atomTypeId(7) == 1);
   TEST_ASSERT(species.atomTypeId(8) == 2);
   TEST_ASSERT(species.atomTypeId(9) == 2);
   TEST_ASSERT(species.atomTypeId(10) == 2);
   TEST_ASSERT(species.atomTypeId(11) == 2);
   if (verbose() > 1) { 
      std::cout << std::endl << "Atom Types:"; 
      for (int i = 0; i < species.nAtom(); ++i) {
         std::cout << std::endl << i << "  "
                   << species.atomTypeId(i);
      }
   }
  
   // Bonds
   TEST_ASSERT(species.nBond() == species.nAtom() - 1);
   const SpeciesGroup<2>* bondPtr;
   if (verbose() > 1) { std::cout << std::endl << "Bonds:"; }
   for (int i = 0; i < species.nBond(); ++i) {
      bondPtr = &species.speciesBond(i);
      TEST_ASSERT(bondPtr->atomId(0) == i);
      TEST_ASSERT(bondPtr->atomId(1) == i+1);
      TEST_ASSERT(bondPtr->typeId() == 0);
      if (verbose() > 1) {
         std::cout << std::endl << i << " "
                   << bondPtr->atomId(0) << " "
                   << bondPtr->atomId(1) << " "
                   << bondPtr->typeId();
      }
   }

   const Species::AtomBondIdArray* bondIdArrayPtr;
   for (int i = 0; i < species.nAtom(); ++i) {
      bondIdArrayPtr = &species.atomBondIds(i);
      if (i == 0) {
         TEST_ASSERT(bondIdArrayPtr->size() == 1);
         TEST_ASSERT((*bondIdArrayPtr)[0] == 0);
      } else
      if (i == species.nAtom() - 1) {
         TEST_ASSERT(bondIdArrayPtr->size() == 1);
         TEST_ASSERT( (*bondIdArrayPtr)[0] == species.nAtom()-2);
      } else {
         TEST_ASSERT(bondIdArrayPtr->size() == 2);
         TEST_ASSERT( (*bondIdArrayPtr)[0] == i-1);
         TEST_ASSERT( (*bondIdArrayPtr)[1] == i);
      }
   }

   #ifdef SIMP_ANGLE
   TEST_ASSERT(species.nAngle() == species.nAtom() - 2);
   const SpeciesGroup<3>* anglePtr;
   if (verbose() > 1) { std::cout << std::endl << "Angles:"; }
   for (int i = 0; i < species.nAngle(); ++i) {
      anglePtr = &species.speciesAngle(i);
      TEST_ASSERT(anglePtr->atomId(0) == i);
      TEST_ASSERT(anglePtr->atomId(1) == i+1);
      TEST_ASSERT(anglePtr->atomId(2) == i+2);
      TEST_ASSERT(anglePtr->typeId() == 0);
      if (verbose() > 1) {
         std::cout << std::endl << i << "  "
                   << anglePtr->atomId(0) << "  "
                   << anglePtr->atomId(1) << "  "
                   << anglePtr->atomId(2) << "  "
                   << anglePtr->typeId();
      }
   }

   const Species::AtomAngleIdArray* angleIdArrayPtr;
   for (int i = 0; i < species.nAtom(); ++i) {
      angleIdArrayPtr = &species.atomAngleIds(i);
      if (i == 0) {
         TEST_ASSERT(angleIdArrayPtr->size() == 1);
         TEST_ASSERT((*angleIdArrayPtr)[0] == i);
      } else
      if (i == 1) {
         TEST_ASSERT(angleIdArrayPtr->size() == 2);
         TEST_ASSERT((*angleIdArrayPtr)[0] == i-1);
         TEST_ASSERT((*angleIdArrayPtr)[1] == i);
      } else
      if (i == species.nAtom() - 1) {
         TEST_ASSERT(angleIdArrayPtr->size() == 1);
         TEST_ASSERT( (*angleIdArrayPtr)[0] == i-2);
      } else
      if (i == species.nAtom() - 2) {
         TEST_ASSERT(angleIdArrayPtr->size() == 2);
         TEST_ASSERT( (*angleIdArrayPtr)[0] == i-2);
         TEST_ASSERT( (*angleIdArrayPtr)[1] == i-1);
      } else {
         TEST_ASSERT(angleIdArrayPtr->size() == 3);
         TEST_ASSERT( (*angleIdArrayPtr)[0] == i-2);
         TEST_ASSERT( (*angleIdArrayPtr)[1] == i-1);
         TEST_ASSERT( (*angleIdArrayPtr)[2] == i);
      }
   }
   #endif
}

TEST_BEGIN(MultiblockTest)
TEST_ADD(MultiblockTest, testConstructor)
TEST_ADD(MultiblockTest, testReadParam)
TEST_END(MultiblockTest)

#endif
