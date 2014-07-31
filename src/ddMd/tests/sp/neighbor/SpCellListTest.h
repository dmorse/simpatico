#ifndef DDMD_SP_CELL_LIST_TEST_H
#define DDMD_SP_CELL_LIST_TEST_H

#include <ddMd/sp/neighbor/SpCellList.h>
#include <ddMd/sp/neighbor/SpCell.h>
#include <ddMd/sp/chemistry/SpAtom.h>
#include <util/containers/DPArray.h>
#include <util/containers/DArray.h>
#include <util/space/Vector.h>
#include <util/space/IntVector.h>
#include <util/random/Random.h>
#include <util/format/Int.h>
#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <sstream>

using namespace Util;
using namespace DdMd;

class SpCellListTest : public UnitTest 
{

private:

   SpCellList spCellList;

public:

   void setUp()
   {}

   void tearDown()
   {}

   void testAllocate()
   {
      printMethod(TEST_FUNC);
      Vector lengths(4.0, 6.0, 8.0);
      Vector lower(2.0, 0.0, 4.0);
      Vector upper(4.0, 3.0, 8.0);
      Vector cutoffs;
      double cutoff = 1.2 / 3.0;
      for (int i=0; i < 3; ++i) {
         lower[i] = lower[i]/lengths[i];
         upper[i] = upper[i]/lengths[i];
         cutoffs[i] = cutoff/lengths[i];
      }
      spCellList.allocate(10, lower, upper, cutoffs);

      TEST_ASSERT(spCellList.atomCapacity() == 10);
      TEST_ASSERT(spCellList.cellCapacity() == 5*7*10);
   }

   void testMakeGrid()
   {
      printMethod(TEST_FUNC);
      Vector lengths(4.0, 6.0, 8.0);
      Vector lower(2.0, 0.0, 4.0);
      Vector upper(4.0, 3.0, 8.0);
      Vector cutoffs;
      double cutoff = 1.2 / 3.0;
      for (int i=0; i < 3; ++i) {
         lower[i] = lower[i]/lengths[i];
         upper[i] = upper[i]/lengths[i];
         cutoffs[i] = cutoff/lengths[i];
      }
      spCellList.allocate(10, lower, upper, cutoffs);
      spCellList.makeGrid(lower, upper, cutoffs);

      TEST_ASSERT(spCellList.grid().dimension(0) == 5);
      TEST_ASSERT(spCellList.grid().dimension(1) == 7);
      TEST_ASSERT(spCellList.grid().dimension(2) == 10);
      TEST_ASSERT(spCellList.grid().size() == 5*7*10);

      // Test length of linked list
      int n = 0;
      const SpCell* spCellPtr = spCellList.begin();
      while (spCellPtr) {
         spCellPtr = spCellPtr->nextCellPtr();
         ++n;
      }
      TEST_ASSERT(n == 5*7*10);

   }

   // this test needs to be recalculated because right now
   // the grid size is 6 which is less than 27 and will cause
   // an error to be thrown if debugging is enabled.
   void testPlaceAtom()
   {
      printMethod(TEST_FUNC);

      const int    maxNAtom = 10;

      // Create a Boundary
      Vector lengths(4.0, 6.0, 8.0);
      Vector lower(2.0, 0.0, 4.0);
      Vector upper(4.0, 3.0, 8.0);
      Vector r(3.1, 1.6, 5.3);  // (0.775, 0.266, 0.6625)
      Vector cutoffs;
      double cutoff = 1.2;
      for (int i=0; i < 3; ++i) {
         lower[i] = lower[i]/lengths[i];
         upper[i] = upper[i]/lengths[i];
         cutoffs[i] = cutoff/lengths[i];
         r[i] = r[i]/lengths[i];
      }
      spCellList.allocate(10, lower, upper, cutoffs);
      spCellList.makeGrid(lower, upper, cutoffs);

      TEST_ASSERT(spCellList.grid().dimension(0) == 1);
      TEST_ASSERT(spCellList.grid().dimension(1) == 2);
      TEST_ASSERT(spCellList.grid().dimension(2) == 3);
      TEST_ASSERT(spCellList.grid().size() == 6);

      // Allocate and initialize an array of Atoms
      DArray<SpAtom>       spAtoms;
      spAtoms.allocate(maxNAtom);

      int cellId  = spCellList.cellIndexFromPosition(r);
      IntVector p = spCellList.grid().position(cellId);
      TEST_ASSERT(p[0] == 0);
      TEST_ASSERT(p[1] == 1);
      TEST_ASSERT(p[2] == 0);
      TEST_ASSERT(cellId == 3);

      int atomId = 3;
      spAtoms[atomId].position = r;
      spCellList.placeAtom(spAtoms[atomId]);
      spCellList.build();
      TEST_ASSERT(spCellList.isValid());
      spCellList.update();

      int size, capacity;
      for (int i = 0; i < spCellList.grid().size(); ++i) {
         size     = spCellList.cell(i).nAtom();
         capacity = spCellList.cell(i).atomCapacity();
         TEST_ASSERT(size == capacity);
         if (i == cellId) {
            TEST_ASSERT(size == 1);
         } else {
            TEST_ASSERT(size == 0);
         }
      }

   }

   // same with this test
   void testPlaceAtoms()
   {
      printMethod(TEST_FUNC);

      // Add nAtom atoms to the same cell
      const int nAtom = 2;
      Vector r[nAtom] = {
                            Vector(3.0, 1.6, 5.0),
                            Vector(2.5, 1.8, 5.2)
                        };
      int cellId[nAtom];

      // Initialize a Boundary and spCellList geometry
      Vector lengths(4.0, 6.0, 8.0);
      Vector lower(2.0, 0.0, 4.0);
      Vector upper(4.0, 3.0, 8.0);
      Vector cutoffs;
      double cutoff = 1.2;
      for (int i=0; i < 3; ++i) {
         lower[i] = lower[i]/lengths[i];
         upper[i] = upper[i]/lengths[i];
         cutoffs[i] = cutoff/lengths[i];
         for (int j=0; j < nAtom; ++j) {
            r[j][i] = r[j][i]/lengths[i];
         }
      }
      spCellList.allocate(10, lower, upper, cutoffs);
      spCellList.makeGrid(lower, upper, cutoffs);

      // Allocate and initialize an array of nAtom atoms
      DArray<SpAtom>       spAtoms;
      spAtoms.allocate(nAtom);

      // Add nAtom atoms to the same cell
      int i;
      spCellList.clear();
      for (i=0; i < nAtom; i++){
         cellId[i] = spCellList.cellIndexFromPosition(r[i]);
         spAtoms[i].position = r[i];
         spCellList.placeAtom(spAtoms[i]);
      }
      spCellList.build();
      spCellList.update();

      TEST_ASSERT( spCellList.nAtom() == nAtom );
      for (i = 0; i < nAtom; ++i){
         TEST_ASSERT(spCellList.cell(cellId[i]).nAtom() == nAtom);
      }

      try {
         spCellList.isValid();
      } catch (Exception e) {
         TEST_ASSERT(0);
      }

   }

   // and this test
   void testAddRandomAtoms()
   {
      printMethod(TEST_FUNC);

      const int nAtom = 100;

      // Setup CellList
      Vector lengths(4.0, 6.0, 8.0);
      Vector lower(2.0, 0.0, 4.0);
      Vector upper(4.0, 3.0, 8.0);
      Vector cutoffs;
      double cutoff = 1.2;
      for (int i = 0; i < 3; ++i) {
         lower[i] = lower[i]/lengths[i];
         upper[i] = upper[i]/lengths[i];
         cutoffs[i] = cutoff/lengths[i];
      }
      spCellList.allocate(nAtom, lower, upper, cutoffs);
      spCellList.makeGrid(lower, upper, cutoffs);

      TEST_ASSERT(spCellList.grid().dimension(0) == 1);
      TEST_ASSERT(spCellList.grid().dimension(1) == 2);
      TEST_ASSERT(spCellList.grid().dimension(2) == 3);
      TEST_ASSERT(eq(spCellList.cellLength(0), 0.5));
      TEST_ASSERT(eq(spCellList.cellLength(1), 0.25));
      TEST_ASSERT(eq(spCellList.cellLength(2), 1.0/6.0));
      TEST_ASSERT(spCellList.grid().size() == 6);

      // Allocate Atom 
      DArray<SpAtom>       spAtoms;
      spAtoms.allocate(nAtom);

      // Create new random number generator
      Random random;
      random.setSeed(1098640);

      // Place spAtoms at random in extended region
      Vector pos;
      spCellList.clear();

      int i;
      for (i = 0; i < nAtom; ++i) {
         spAtoms[i].typeId = 1;
         for (int j = 0; j < Dimension; ++j) {
            pos[j] = random.uniform(lower[j], upper[j]);
            TEST_ASSERT(pos[j] >= lower[j]);
            TEST_ASSERT(pos[j] <= upper[j]);
         }
         spAtoms[i].position = pos;

         try {
            spCellList.placeAtom(spAtoms[i]);
         }
         catch (Exception e) {
            e.write(std::cout);
            TEST_ASSERT(0);
         }
      }

      spCellList.build();
      spCellList.update();

      int na = 0;
      int nc = 0;
      const SpCell* spCellPtr = spCellList.begin();
      while (spCellPtr) {
         ++nc;
         na += spCellPtr->nAtom();
         spCellPtr = spCellPtr->nextCellPtr();
      }
      TEST_ASSERT(na == nAtom);

      try {
         spCellList.isValid();
      } catch (Exception e) {
         TEST_ASSERT(0);
      }

   }

};

TEST_BEGIN(SpCellListTest)
TEST_ADD(SpCellListTest, testAllocate)
TEST_ADD(SpCellListTest, testMakeGrid)
//TEST_ADD(SpCellListTest, testPlaceAtom)
//TEST_ADD(SpCellListTest, testPlaceAtoms)
//TEST_ADD(SpCellListTest, testAddRandomAtoms)
TEST_END(SpCellListTest)

#endif //ifndef DDMD_SP_CELL_LIST_TEST_H
