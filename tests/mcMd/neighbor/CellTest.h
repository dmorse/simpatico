#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <mcMd/neighbor/Cell.h>
#include <mcMd/neighbor/CellTag.h>
#include <mcMd/chemistry/Atom.h>
#include <util/containers/RArray.h>
#include <util/containers/DArray.h>

using namespace Util;
using namespace McMd;

class CellTest : public UnitTest 
{

private:

   static const int nAtom  = 6;
   static const int cellId = 37;

   Cell            cell;
   DArray<CellTag> cellTags;
   RArray<Atom>    atoms;

public:

   void setUp()
   {
      cellTags.allocate(nAtom);
      Atom::allocate(nAtom, atoms);
   }

   void tearDown()
   {  Atom::deallocate(); }

   void testAddAtom()
   {
      printMethod(TEST_FUNC);

      for (int i=0; i< nAtom; i++) {
         cell.addAtom(cellTags[i], atoms[i], cellId);
         TEST_ASSERT(cellTags[i].cellId  == cellId);
         TEST_ASSERT(cellTags[i].cellPos == i);
         TEST_ASSERT(cell.nAtomCell_ == i+1);
         TEST_ASSERT(cell.firstEmptyPos_ == i+1);
         TEST_ASSERT(cell.firstClearPos_ == i+1);
      }
      TEST_ASSERT(cell.nAtomCell_ == nAtom);
      TEST_ASSERT(cell.firstEmptyPos_ == nAtom);
      TEST_ASSERT(cell.firstClearPos_ == nAtom);

      cell.isValid(cellTags, nAtom, cellId);

   }

   void testAddDeleteAtom()
   {
      printMethod(TEST_FUNC);
      int i;

      for (i=0; i< nAtom; i++) {
         cell.addAtom(cellTags[i], atoms[i], cellId);
      }
      cell.deleteAtom(cellTags[2]);
      TEST_ASSERT(cell.nAtomCell_ == nAtom-1);
      TEST_ASSERT(cell.firstEmptyPos_ == 2);
      TEST_ASSERT(cell.firstClearPos_ == nAtom);
      cell.isValid(cellTags, nAtom, cellId);

      cell.deleteAtom(cellTags[4]);

      TEST_ASSERT(cell.nAtomCell_ == nAtom-2);
      TEST_ASSERT(cell.firstEmptyPos_ == 2);
      TEST_ASSERT(cell.firstClearPos_ == nAtom);
      cell.isValid(cellTags, nAtom, cellId);

      cell.deleteAtom(cellTags[5]);

      TEST_ASSERT(cell.nAtomCell_ == nAtom-3);
      TEST_ASSERT(cell.firstEmptyPos_ == 2);
      TEST_ASSERT(cell.firstClearPos_ == 4);

      cell.isValid(cellTags, nAtom, cellId);

      cell.addAtom(cellTags[2], atoms[2], cellId);

      TEST_ASSERT(cell.nAtomCell_ == nAtom-2);
      TEST_ASSERT(cell.firstEmptyPos_ == 4);
      TEST_ASSERT(cell.firstClearPos_ == 4);
      TEST_ASSERT(cellTags[2].cellId == cellId);
      TEST_ASSERT(cellTags[2].cellPos == 2);
      cell.isValid(cellTags, nAtom, cellId);

   }

};

TEST_BEGIN(CellTest)
TEST_ADD(CellTest, testAddAtom)
TEST_ADD(CellTest, testAddDeleteAtom)
TEST_END(CellTest)


