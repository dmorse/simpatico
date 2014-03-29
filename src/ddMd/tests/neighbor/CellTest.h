#include <ddMd/neighbor/Cell.h>
#include <ddMd/neighbor/CellAtom.h>
#include <ddMd/chemistry/Atom.h>
#include <ddMd/chemistry/AtomArray.h>
#include <util/containers/RArray.h>
#include <util/containers/DArray.h>

#ifdef UTIL_MPI
#ifndef TEST_MPI
#define TEST_MPI
#endif
#endif

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

using namespace Util;
using namespace DdMd;

class CellTest : public UnitTest 
{

private:

   static const int nAtom  = 10;

   Cell  cell;
   AtomArray  atoms;
   DArray<CellAtom>  cellAtoms;

public:

   void setUp()
   {
      atoms.allocate(nAtom);
      cellAtoms.allocate(nAtom);
   }

   void tearDown()
   {}  

   void testInitialize() 
   {
      printMethod(TEST_FUNC);

      TEST_ASSERT(cell.nAtom() == 0);
      TEST_ASSERT(cell.atomCapacity() == 0);
      for (int i = 0; i < nAtom; ++i) {
         cell.incrementCapacity();
      }
      cell.initialize(&(cellAtoms[0]));
      TEST_ASSERT(cell.nAtom() == 0);
      TEST_ASSERT(cell.atomCapacity() == nAtom);
   }

   void testAddAtoms() 
   {
      printMethod(TEST_FUNC);

      cell.clear();
      for (int i = 0; i < nAtom; ++i) {
         cell.incrementCapacity();
      }
      cell.initialize(&(cellAtoms[0]));
      TEST_ASSERT(cell.nAtom() == 0);
      TEST_ASSERT(cell.atomCapacity() == nAtom);
      for (int i = 0; i < nAtom; ++i) {
         TEST_ASSERT(cell.nAtom() == i);
         cell.append(&(atoms[i]));
         TEST_ASSERT(cell.atomPtr(i) == &(cellAtoms[i]));
      }
      TEST_ASSERT(cell.nAtom() == nAtom);
      TEST_ASSERT(cell.atomCapacity() == nAtom);
      for (int i = 0; i < nAtom; ++i) {
         TEST_ASSERT(cell.atomPtr(i) == &(cellAtoms[i]));
         TEST_ASSERT(cellAtoms[i].ptr() == &(atoms[i]));
      }
      
   }

};

TEST_BEGIN(CellTest)
TEST_ADD(CellTest, testInitialize)
TEST_ADD(CellTest, testAddAtoms)
TEST_END(CellTest)
