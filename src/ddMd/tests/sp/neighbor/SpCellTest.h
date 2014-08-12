#ifndef DDMD_SP_CELL_TEST_H
#define DDMD_SP_CELL_TEST_H

#include <ddMd/sp/neighbor/SpCell.h>
#include <ddMd/sp/chemistry/SpAtom.h>
#include <util/containers/DArray.h>
#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

using namespace Util;
using namespace DdMd;

class SpCellTest : public UnitTest
{
private:

   static const int nAtom  = 10;

   SpCell               spCell;
   DArray<SpCellAtom>   spCellAtoms;
   DArray<SpAtom>       spAtoms;

public:

   void setUp()
   {
      spAtoms.allocate(nAtom);
      spCellAtoms.allocate(nAtom);
   }

   void tearDown()
   {}

   void testInitialize() 
   {
      printMethod(TEST_FUNC);

      TEST_ASSERT(spCell.nAtom() == 0);
      TEST_ASSERT(spCell.atomCapacity() == 0);
      for (int i = 0; i < nAtom; ++i) {
         spCell.incrementCapacity();
      }
      spCell.initialize(&(spCellAtoms[0]));
      TEST_ASSERT(spCell.nAtom() == 0);
      TEST_ASSERT(spCell.atomCapacity() == nAtom);
   }

   void testAddAtoms() 
   {
      printMethod(TEST_FUNC);

      spCell.clear();
      for (int i = 0; i < nAtom; ++i) {
         spCell.incrementCapacity();
      }
      spCell.initialize(&(spCellAtoms[0]));
      TEST_ASSERT(spCell.nAtom() == 0);
      TEST_ASSERT(spCell.atomCapacity() == nAtom);
      for (int i = 0; i < nAtom; ++i) {
         TEST_ASSERT(spCell.nAtom() == i);
         spCell.append(&(spAtoms[i]));
         TEST_ASSERT(spCell.atomPtr(i) == &(spCellAtoms[i]));
      }
      TEST_ASSERT(spCell.nAtom() == nAtom);
      TEST_ASSERT(spCell.atomCapacity() == nAtom);
      for (int i = 0; i < nAtom; ++i) {
         TEST_ASSERT(spCell.atomPtr(i) == &(spCellAtoms[i]));
         TEST_ASSERT(spCellAtoms[i].ptr() == &(spAtoms[i]));
      }

   }
};

TEST_BEGIN(SpCellTest)
TEST_ADD(SpCellTest, testInitialize)
TEST_ADD(SpCellTest, testAddAtoms)
TEST_END(SpCellTest)

#endif //ifndef DDMD_SP_CELL_TEST_H
