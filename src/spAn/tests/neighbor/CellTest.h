#ifndef SPAN_CELL_TEST_H
#define SPAN_CELL_TEST_H

#include <spAn/neighbor/Cell.h>
#include <spAn/chemistry/Atom.h>
#include <util/containers/DArray.h>
#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

using namespace Util;
using namespace SpAn;

class CellTest : public UnitTest
{
private:

   static const int nAtom  = 10;

   Cell               spCell;
   DArray<CellAtom>   spCellAtoms;
   DArray<Atom>       spAtoms;

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

TEST_BEGIN(CellTest)
TEST_ADD(CellTest, testInitialize)
TEST_ADD(CellTest, testAddAtoms)
TEST_END(CellTest)

#endif //ifndef SPAN_CELL_TEST_H
