#ifndef SPAN_CELL_LIST_TEST_H
#define SPAN_CELL_LIST_TEST_H

#include <spAn/neighbor/CellList.h>
#include <spAn/neighbor/Cell.h>
#include <spAn/chemistry/Atom.h>
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
using namespace SpAn;

class CellListTest : public UnitTest 
{

private:

   CellList spCellList;

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
      const Cell* spCellPtr = spCellList.begin();
      while (spCellPtr) {
         spCellPtr = spCellPtr->nextCellPtr();
         ++n;
      }
      TEST_ASSERT(n == 5*7*10);

   }

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
      double cutoff = 1.2 / 3.0;
      for (int i=0; i < 3; ++i) {
         lower[i] = lower[i]/lengths[i];
         upper[i] = upper[i]/lengths[i];
         cutoffs[i] = cutoff/lengths[i];
         r[i] = r[i]/lengths[i];
      }
      spCellList.allocate(10, lower, upper, cutoffs);
      spCellList.makeGrid(lower, upper, cutoffs);

      TEST_ASSERT(spCellList.grid().dimension(0) == 5);
      TEST_ASSERT(spCellList.grid().dimension(1) == 7);
      TEST_ASSERT(spCellList.grid().dimension(2) == 10);
      TEST_ASSERT(spCellList.grid().size() == 5*7*10);

      // Allocate and initialize an array of Atoms
      DArray<Atom>       spAtoms;
      spAtoms.allocate(maxNAtom);

      int cellId  = spCellList.cellIndexFromPosition(r);
      IntVector p = spCellList.grid().position(cellId);
      TEST_ASSERT(p[0] == 2);
      TEST_ASSERT(p[1] == 3);
      TEST_ASSERT(p[2] == 3);
      TEST_ASSERT(cellId == 173);

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
      double cutoff = 1.2 / 3.0;
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
      DArray<Atom>       spAtoms;
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
         TEST_ASSERT(spCellList.cell(cellId[i]).nAtom() == 1);
      }

      try {
         spCellList.isValid();
      } catch (Exception e) {
         TEST_ASSERT(0);
      }

   }

   void testAddRandomAtoms()
   {
      printMethod(TEST_FUNC);

      const int nAtom = 100;

      // Setup CellList
      Vector lengths(4.0, 6.0, 8.0);
      Vector lower(2.0, 0.0, 4.0);
      Vector upper(4.0, 3.0, 8.0);
      Vector cutoffs;
      double cutoff = 1.2 / 3.0;
      for (int i = 0; i < 3; ++i) {
         lower[i] = lower[i]/lengths[i];
         upper[i] = upper[i]/lengths[i];
         cutoffs[i] = cutoff/lengths[i];
      }
      spCellList.allocate(nAtom, lower, upper, cutoffs);
      spCellList.makeGrid(lower, upper, cutoffs);

      TEST_ASSERT(spCellList.grid().dimension(0) == 5);
      TEST_ASSERT(spCellList.grid().dimension(1) == 7);
      TEST_ASSERT(spCellList.grid().dimension(2) == 10);
      TEST_ASSERT(eq(spCellList.cellLength(0), 0.5 / 5.0));
      TEST_ASSERT(eq(spCellList.cellLength(1), 0.5 / 7.0));
      TEST_ASSERT(eq(spCellList.cellLength(2), 0.5 / 10.0));
      TEST_ASSERT(spCellList.grid().size() == 5*7*10);

      // Allocate Atom 
      DArray<Atom>       spAtoms;
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
      const Cell* spCellPtr = spCellList.begin();
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

   void testGetNeighbors()
   {
      printMethod(TEST_FUNC);

      // Define cutoff
      double cutoff = 1.20 / 3.0;
      double pairCutoffSq = cutoff - 1.0E-10;
      pairCutoffSq = pairCutoffSq*pairCutoffSq;

      // Setup domain
      Vector lengths(4.0, 6.0, 8.0);
      Vector lower(2.0, 0.0, 4.0);
      Vector upper(4.0, 3.0, 8.0);
      Vector cutoffs;
      for (int i=0; i < Dimension; ++i) {
         lower[i] = lower[i]/lengths[i];
         upper[i] = upper[i]/lengths[i];
         cutoffs[i] = cutoff/lengths[i];
      }

      // Setup cell list grid
      const int nAtom = 200;
      spCellList.allocate(nAtom, lower, upper, cutoffs);
      spCellList.makeGrid(lower, upper, cutoffs);

      TEST_ASSERT(spCellList.grid().dimension(0) == 5);
      TEST_ASSERT(spCellList.grid().dimension(1) == 7);
      TEST_ASSERT(spCellList.grid().dimension(2) == 10);
      TEST_ASSERT(eq(spCellList.cellLength(0), 0.5/5.0));
      TEST_ASSERT(eq(spCellList.cellLength(1), 0.5/7.0));
      TEST_ASSERT(eq(spCellList.cellLength(2), 0.5/10.0));
      TEST_ASSERT(spCellList.grid().size() == 5*7*10);

      // Allocate Atom 
      DArray<Atom> spAtoms;
      DPArray<Atom> locals;
      spAtoms.allocate(nAtom);
      locals.allocate(nAtom);

      // Create new random number generator
      Random random;
      random.setSeed(1098640);

      // Place atoms at random in extended region
      Vector pos;
      int nLocal = 0;
      spCellList.clear();
      for (int i = 0; i < nAtom; ++i) {
         spAtoms[i].typeId = 1;
         for (int j = 0; j < Dimension; ++j) {
            pos[j] = random.uniform(lower[j], upper[j]);
            TEST_ASSERT(pos[j] >= lower[j]);
            TEST_ASSERT(pos[j] <= upper[j]);
         }
         spAtoms[i].position = pos;
         ++nLocal;
         locals.append(spAtoms[i]);

         TEST_ASSERT(nLocal == locals.size());

         try {
            spCellList.placeAtom(spAtoms[i]);
         }
         catch (Exception e) {
            if (isIoProcessor()) {
               e.write(std::cout);
            }
            TEST_ASSERT(0);
         }

      }

      spCellList.build();

      try {
         spCellList.isValid();
      }
      catch (Exception e) {
         TEST_ASSERT(0);
      }

      // Check that # atoms in local cells = # of local atoms
      int na = 0;
      int nc = 0;
      const Cell* spCellPtr = spCellList.begin();
      while (spCellPtr) {
         ++nc;
         na += spCellPtr->nAtom();
         spCellPtr = spCellPtr->nextCellPtr();
      }
      TEST_ASSERT(na == locals.size());

      // Transform coordinates to Cartesian
      for (int i = 0; i < nAtom; ++i) {
         for (int j = 0; j < Dimension; ++j) {
            spAtoms[i].position[j] *= lengths[j];
         }
      }
      spCellList.update();

      // Find all neighbor pairs within a cutoff (use cell list)
      Cell::NeighborArray neighbors;
      CellAtom* spCellAtomPtr1;
      CellAtom* spCellAtomPtr2;
      Vector dr;
      int    nn;      // number of neighbors in a cell
      int    np = 0;  // Number of pairs within cutoff
      spCellPtr = spCellList.begin();
      while (spCellPtr) {
         spCellPtr->getNeighbors(neighbors);
         na = spCellPtr->nAtom();
         nn = neighbors.size();
         for (int i = 0; i < na; ++i) {
            spCellAtomPtr1 = neighbors[i];
            for (int j = 0; j < na; ++j) {
               spCellAtomPtr2 = neighbors[j];
               if (spCellAtomPtr2 > spCellAtomPtr1) {
                  dr.subtract(spCellAtomPtr2->position(), spCellAtomPtr1->position());
                  if (dr.square() < pairCutoffSq) {
                     ++np;
                  }
               }
            }
            for (int j = na; j < nn; ++j) {
               spCellAtomPtr2 = neighbors[j];
               dr.subtract(spCellAtomPtr2->position(), spCellAtomPtr1->position());
               if (dr.square() < pairCutoffSq) {
                  ++np;
               }
            }
         }
         spCellPtr = spCellPtr->nextCellPtr();
      }

      // Count neighbor pairs directly (N^2 loop)
      Atom* spAtomPtr1;
      Atom* spAtomPtr2;
      int nq = 0;
      for (int i = 0; i < locals.size(); ++i) {
         spAtomPtr1 = &locals[i];
         for (int j = 0; j < locals.size(); ++j) {
            spAtomPtr2 = &locals[j];
            if (spAtomPtr2 > spAtomPtr1) {
               dr.subtract(spAtomPtr2->position, spAtomPtr1->position);
               if (dr.square() < pairCutoffSq) {
                  ++nq;
               }
            }
         }
      }

      // Assert that number of neighobr pairs is same by either method.
      TEST_ASSERT(np == nq);
   }

};

TEST_BEGIN(CellListTest)
TEST_ADD(CellListTest, testAllocate)
TEST_ADD(CellListTest, testMakeGrid)
TEST_ADD(CellListTest, testPlaceAtom)
TEST_ADD(CellListTest, testPlaceAtoms)
TEST_ADD(CellListTest, testAddRandomAtoms)
TEST_ADD(CellListTest, testGetNeighbors)
TEST_END(CellListTest)

#endif //ifndef SPAN_CELL_LIST_TEST_H
