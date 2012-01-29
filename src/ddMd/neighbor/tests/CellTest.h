#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <ddMd/neighbor/Cell.h>
#include <ddMd/chemistry/Atom.h>
#include <ddMd/chemistry/AtomArray.h>
#include <util/containers/RArray.h>
#include <util/containers/DArray.h>

using namespace Util;
using namespace DdMd;

class CellTest : public UnitTest 
{

private:

   typedef Atom* Handle;

   static const int nAtom  = 10;

   Cell         cell;
   AtomArray      atoms;
   DArray<Atom*>  handles;

public:

   void setUp()
   {
      atoms.allocate(nAtom);
      handles.allocate(nAtom);
      for (int i = 0; i < nAtom; i++) {
         handles[i] = &atoms[i];
      }
   }

   void tearDown()
   {}  

   void testInitialize() 
   {
      printMethod(TEST_FUNC);

      TEST_ASSERT(cell.nAtom() == 0);
      TEST_ASSERT(cell.atomCapacity() == 0);
      for (int i = 0; i < nAtom; i++) {
         cell.incrementCapacity();
      }
      cell.initialize(&handles[0]);
      TEST_ASSERT(cell.nAtom() == 0);
      TEST_ASSERT(cell.atomCapacity() == nAtom);
   }

   void testAddAtoms() 
   {
      printMethod(TEST_FUNC);

      for (int i = 0; i < nAtom; i++) {
         cell.incrementCapacity();
      }
      cell.initialize(&handles[0]);
      TEST_ASSERT(cell.nAtom() == 0);
      TEST_ASSERT(cell.atomCapacity() == nAtom);
      for (int i = 0; i < nAtom; ++i) {
         cell.append(handles[i]);
      }
      TEST_ASSERT(cell.nAtom() == nAtom);
      TEST_ASSERT(cell.atomCapacity() == nAtom);
      for (int i = 0; i < nAtom; ++i) {
         TEST_ASSERT(cell.atomPtr(i) == handles[i]);
         TEST_ASSERT(cell.atomPtr(i) == &atoms[i]);
      }
      
   }

};

TEST_BEGIN(CellTest)
TEST_ADD(CellTest, testInitialize)
TEST_ADD(CellTest, testAddAtoms)
TEST_END(CellTest)


