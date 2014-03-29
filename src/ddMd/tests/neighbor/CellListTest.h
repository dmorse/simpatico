#ifndef DDMD_CELL_LIST_TEST_H
#define DDMD_CELL_LIST_TEST_H

#include <ddMd/neighbor/CellList.h>
#include <ddMd/chemistry/Atom.h>
#include <ddMd/chemistry/AtomArray.h>
#include <util/containers/DPArray.h>
#include <util/space/Vector.h>
#include <util/random/Random.h>
#include <util/format/Int.h>

#ifdef UTIL_MPI
#ifndef TEST_MPI
#define TEST_MPI
#endif
#endif

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <iostream>

using namespace Util;
using namespace DdMd;

class CellListTest : public UnitTest 
{

private:

   CellList cellList;

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
      double cutoff = 1.2;
      for (int i=0; i < 3; ++i) {
         lower[i] = lower[i]/lengths[i]; 
         upper[i] = upper[i]/lengths[i]; 
         cutoffs[i] = cutoff/lengths[i]; 
      }
      cellList.allocate(10, lower, upper, cutoffs);

      TEST_ASSERT(cellList.atomCapacity() == 10);
      TEST_ASSERT(cellList.cellCapacity() == 60);
   }

   void testAllocate2()
   {
      printMethod(TEST_FUNC);
      Vector lengths(4.0, 6.0, 8.0);
      Vector lower(2.0, 0.0, 4.0);
      Vector upper(4.0, 3.0, 8.0);
      Vector cutoffs;
      double cutoff = 1.19;
      for (int i=0; i < 3; ++i) {
         lower[i] = lower[i]/lengths[i]; 
         upper[i] = upper[i]/lengths[i]; 
         cutoffs[i] = cutoff/lengths[i]; 
      }
      int nCellCut = 3;
      cellList.allocate(10, lower, upper, cutoffs, nCellCut);

      TEST_ASSERT(cellList.atomCapacity() == 10);
      TEST_ASSERT(cellList.cellCapacity() == 11*13*16);
   }

   void testMakeGrid()
   {
      printMethod(TEST_FUNC);
      Vector lengths(4.0, 6.0, 8.0);
      Vector lower(2.0, 0.0, 4.0);
      Vector upper(4.0, 3.0, 8.0);
      Vector cutoffs;
      double cutoff = 1.2;
      for (int i=0; i < 3; ++i) {
         lower[i] = lower[i]/lengths[i]; 
         upper[i] = upper[i]/lengths[i]; 
         cutoffs[i] = cutoff/lengths[i]; 
      }
      cellList.allocate(10, lower, upper, cutoffs);
      cellList.makeGrid(lower, upper, cutoffs);

      TEST_ASSERT(cellList.grid().dimension(0) == 3);
      TEST_ASSERT(cellList.grid().dimension(1) == 4);
      TEST_ASSERT(cellList.grid().dimension(2) == 5);
      TEST_ASSERT(cellList.grid().size() == 60);

      // Test length of linked list
      int n = 0;
      const Cell* cellPtr = cellList.begin();
      while (cellPtr) {
         cellPtr = cellPtr->nextCellPtr();
         ++n;
      }
      TEST_ASSERT(n == 6);

   }

   void testMakeGrid2()
   {
      printMethod(TEST_FUNC);
      Vector lengths(4.0, 6.0, 8.0);
      Vector lower(2.0, 0.0, 4.0);
      Vector upper(4.0, 3.0, 8.0);
      Vector cutoffs;
      double cutoff = 1.19;
      for (int i=0; i < 3; ++i) {
         lower[i] = lower[i]/lengths[i]; 
         upper[i] = upper[i]/lengths[i]; 
         cutoffs[i] = cutoff/lengths[i]; 
      }
      int nCellCut = 3;
      cellList.allocate(10, lower, upper, cutoffs, nCellCut);
      cellList.makeGrid(lower, upper, cutoffs, nCellCut);

      TEST_ASSERT(cellList.grid().dimension(0) == 11);
      TEST_ASSERT(cellList.grid().dimension(1) == 13);
      TEST_ASSERT(cellList.grid().dimension(2) == 16);
      TEST_ASSERT(cellList.grid().size() == 11*13*16);


      // Loop over local cells, test length of linked list.
      const Cell* cellPtr = cellList.begin();
      int n = 0;
      while (cellPtr) {
         cellPtr = cellPtr->nextCellPtr();
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
      Vector r(3.1, 1.6, 5.3);
      Vector cutoffs;
      double cutoff = 1.2;
      for (int i=0; i < 3; ++i) {
         lower[i] = lower[i]/lengths[i]; 
         upper[i] = upper[i]/lengths[i]; 
         cutoffs[i] = cutoff/lengths[i]; 
         r[i] = r[i]/lengths[i];
      }
      cellList.allocate(10, lower, upper, cutoffs);
      cellList.makeGrid(lower, upper, cutoffs);

      TEST_ASSERT(cellList.grid().dimension(0) == 3);
      TEST_ASSERT(cellList.grid().dimension(1) == 4);
      TEST_ASSERT(cellList.grid().dimension(2) == 5);
      TEST_ASSERT(cellList.grid().size() == 60);

      // Allocate and initialize an array of Atoms
      AtomArray atoms;
      atoms.allocate(maxNAtom);

      int cellId  = cellList.cellIndexFromPosition(r);
      IntVector p = cellList.grid().position(cellId);
      TEST_ASSERT(p[0] == 1);
      TEST_ASSERT(p[1] == 2);
      TEST_ASSERT(p[1] == 2);
      TEST_ASSERT(cellId == 31);

      int atomId = 3;
      atoms[atomId].position() = r;
      cellList.placeAtom(atoms[atomId]);
      cellList.build();
      TEST_ASSERT(cellList.isValid());
      cellList.update();

      int size, capacity;
      for (int i = 0; i < cellList.grid().size(); ++i) {
         size     = cellList.cell(i).nAtom();
         capacity = cellList.cell(i).atomCapacity();
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

      // Initialize a Boundary and CellList geometry
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
      cellList.allocate(10, lower, upper, cutoffs);
      cellList.makeGrid(lower, upper, cutoffs);

      // Allocate and initialize an array of nAtom atoms
      AtomArray atoms;
      atoms.allocate(nAtom);

      // Add nAtom atoms to the same cell
      int i;
      cellList.clear();
      for (i=0; i < nAtom; i++){
         cellId[i] = cellList.cellIndexFromPosition(r[i]);
         atoms[i].position() = r[i];
         cellList.placeAtom(atoms[i]);
      }
      cellList.build();
      cellList.update();

      TEST_ASSERT( cellList.nAtom() == nAtom );
      for (i = 0; i < nAtom; ++i){
         TEST_ASSERT(cellList.cell(cellId[i]).nAtom() == nAtom);
      }

      try { 
         cellList.isValid(); 
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
      double cutoff = 1.2;
      for (int i=0; i < 3; ++i) {
         lower[i] = lower[i]/lengths[i]; 
         upper[i] = upper[i]/lengths[i]; 
         cutoffs[i] = cutoff/lengths[i]; 
      }
      cellList.allocate(nAtom, lower, upper, cutoffs);
      cellList.makeGrid(lower, upper, cutoffs);

      TEST_ASSERT(cellList.grid().dimension(0) == 3);
      TEST_ASSERT(cellList.grid().dimension(1) == 4);
      TEST_ASSERT(cellList.grid().dimension(2) == 5);
      TEST_ASSERT(eq(cellList.cellLength(0), 0.5));
      TEST_ASSERT(eq(cellList.cellLength(1), 0.25));
      TEST_ASSERT(eq(cellList.cellLength(2), 1.0/6.0));
      TEST_ASSERT(cellList.grid().size() == 60);

      Vector lowerGhost;
      Vector upperGhost;
      int    i, j;
      for (i = 0; i < Dimension; ++i) {
          lowerGhost[i] = lower[i] - cutoffs[i];
          upperGhost[i] = upper[i] + cutoffs[i];
      }

      // Allocate Atom 
      AtomArray atoms;
      atoms.allocate(nAtom);

      // Create new random number generator
      Random random;
      random.setSeed(1098640);

      // Place atoms at random in extended region
      Vector pos;
      int ic;
      int nLocal = 0;
      int nGhost = 0;
      bool ghost;
      cellList.clear();
      for (i = 0; i < nAtom; ++i) {
         atoms[i].setTypeId(1);
         ghost = false;
         for (j = 0; j < Dimension; ++j) {
            pos[j] = random.uniform(lowerGhost[j], upperGhost[j]);
            if (pos[j] < lower[j]) 
               ghost = true;
            if (pos[j] > upper[j]) 
               ghost = true;
            TEST_ASSERT(pos[j] >= lowerGhost[j]);
            TEST_ASSERT(pos[j] <= upperGhost[j]);
         }
         atoms[i].position() = pos;
         if (ghost) {
            ++nGhost;
         } else {
            ++nLocal;
         }
         ic = cellList.cellIndexFromPosition(pos);
         TEST_ASSERT( ghost == cellList.cell(ic).isGhostCell() );

         try {
            cellList.placeAtom(atoms[i]);
         } 
         catch (Exception e) {
            e.write(std::cout);
            TEST_ASSERT(0);
         }

      }
      cellList.build();
      cellList.update();

      int na = 0;
      int nc = 0;
      const Cell* cellPtr = cellList.begin();
      while (cellPtr) {
         ++nc;
         na += cellPtr->nAtom();
         cellPtr = cellPtr->nextCellPtr();
      }
      TEST_ASSERT(na == nLocal);

      try { 
         cellList.isValid(); 
      } catch (Exception e) {
         TEST_ASSERT(0);
      }

   }

   void testAddRandomAtoms2()
   {
      printMethod(TEST_FUNC);

      const int nAtom = 100;

      // Setup CellList
      Vector lengths(4.0, 6.0, 8.0);
      Vector lower(2.0, 0.0, 4.0);
      Vector upper(4.0, 3.0, 8.0);
      Vector cutoffs;
      double cutoff = 1.19;
      for (int i=0; i < 3; ++i) {
         lower[i] = lower[i]/lengths[i]; 
         upper[i] = upper[i]/lengths[i]; 
         cutoffs[i] = cutoff/lengths[i]; 
      }
      int nCellCut = 3;
      cellList.allocate(nAtom, lower, upper, cutoffs, nCellCut);
      cellList.makeGrid(lower, upper, cutoffs, nCellCut);

      TEST_ASSERT(cellList.grid().dimension(0) == 11);
      TEST_ASSERT(cellList.grid().dimension(1) == 13);
      TEST_ASSERT(cellList.grid().dimension(2) == 16);
      TEST_ASSERT(eq(cellList.cellLength(0), 0.5/5.0));
      TEST_ASSERT(eq(cellList.cellLength(1), 0.5/7.0));
      TEST_ASSERT(eq(cellList.cellLength(2), 0.5/10.0));
      TEST_ASSERT(cellList.grid().size() == 11*13*16);

      Vector lowerGhost;
      Vector upperGhost;
      int    i, j;
      for (i = 0; i < Dimension; ++i) {
         lowerGhost[i] = lower[i] - cutoffs[i];
         upperGhost[i] = upper[i] + cutoffs[i];
      }

      // Allocate Atom 
      AtomArray atoms;
      atoms.allocate(nAtom);

      // Create new random number generator
      Random random;
      random.setSeed(1098640);

      // Place atoms at random in extended region
      Vector pos;
      int ic;
      int nLocal = 0;
      int nGhost = 0;
      bool ghost;
      cellList.clear();
      for (i = 0; i < nAtom; ++i) {
         atoms[i].setTypeId(1);
         ghost = false;
         for (j = 0; j < Dimension; ++j) {
            pos[j] = random.uniform(lowerGhost[j], upperGhost[j]);
            if (pos[j] < lower[j]) 
               ghost = true;
            if (pos[j] > upper[j]) 
               ghost = true;
            TEST_ASSERT(pos[j] >= lowerGhost[j]);
            TEST_ASSERT(pos[j] <= upperGhost[j]);
         }
         atoms[i].position() = pos;
         if (ghost) {
            ++nGhost;
         } else {
            ++nLocal;
         }
         ic = cellList.cellIndexFromPosition(pos);
         TEST_ASSERT(ghost == cellList.cell(ic).isGhostCell());

         try {
            cellList.placeAtom(atoms[i]);
         } 
         catch (Exception e) {
            e.write(std::cout);
            TEST_ASSERT(0);
         }

      }
      cellList.build();
      cellList.update();

      int na = 0;
      int nc = 0;
      const Cell* cellPtr = cellList.begin();
      while (cellPtr) {
         ++nc;
         na += cellPtr->nAtom();
         cellPtr = cellPtr->nextCellPtr();
      }
      TEST_ASSERT(na == nLocal);

      try { 
         cellList.isValid(); 
      } catch (Exception e) {
         TEST_ASSERT(0);
      }

   }

   void testGetNeighbors()
   {
      printMethod(TEST_FUNC);

      // Define cutoff
      double cutoff = 1.19;
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
      cellList.allocate(nAtom, lower, upper, cutoffs);
      cellList.makeGrid(lower, upper, cutoffs);

      TEST_ASSERT(cellList.grid().dimension(0) == 3);
      TEST_ASSERT(cellList.grid().dimension(1) == 4);
      TEST_ASSERT(cellList.grid().dimension(2) == 5);
      TEST_ASSERT(eq(cellList.cellLength(0), 1.0/2.0));
      TEST_ASSERT(eq(cellList.cellLength(1), 1.0/4.0));
      TEST_ASSERT(eq(cellList.cellLength(2), 1.0/6.0));
      TEST_ASSERT(cellList.grid().size() == 60);

      Vector lowerGhost;
      Vector upperGhost;
      int i, j;
      for (i = 0; i < Dimension; ++i) {
         lowerGhost[i] = lower[i] - cutoffs[i];
         upperGhost[i] = upper[i] + cutoffs[i];
      }

      // Allocate Atom 
      AtomArray atoms;
      DPArray<Atom> locals;
      DPArray<Atom> ghosts;
      atoms.allocate(nAtom);
      locals.allocate(nAtom);
      ghosts.allocate(nAtom);

      // Create new random number generator
      Random random;
      random.setSeed(1098640);

      // Place atoms at random in extended region
      Vector pos;
      int ic;
      int nLocal = 0;
      int nGhost = 0;
      bool ghost;
      cellList.clear();
      for (i = 0; i < nAtom; ++i) {
         atoms[i].setTypeId(1);
         ghost = false;
         for (j = 0; j < Dimension; ++j) {
            pos[j] = random.uniform(lowerGhost[j], upperGhost[j]);
            if (pos[j] < lower[j]) 
               ghost = true;
            if (pos[j] > upper[j]) 
               ghost = true;
            TEST_ASSERT(pos[j] >= lowerGhost[j]);
            TEST_ASSERT(pos[j] <= upperGhost[j]);
         }
         atoms[i].position() = pos;
         if (ghost) {
            ++nGhost;
            ghosts.append(atoms[i]);
         } else {
            ++nLocal;
            locals.append(atoms[i]);
         }
         ic = cellList.cellIndexFromPosition(pos);
         TEST_ASSERT(ghost == cellList.cell(ic).isGhostCell());
         TEST_ASSERT(nGhost == ghosts.size());
         TEST_ASSERT(nLocal == locals.size());

         try {
            cellList.placeAtom(atoms[i]);
         } 
         catch (Exception e) {
            if (isIoProcessor()) {
               e.write(std::cout);
            }
            TEST_ASSERT(0);
         }

      }
      cellList.build();

      try { 
         cellList.isValid(); 
      } catch (Exception e) {
         TEST_ASSERT(0);
      }

      // Check that # atoms in local cells = # of local atoms
      int na = 0;
      int nc = 0;
      const Cell* cellPtr = cellList.begin();
      while (cellPtr) {
         ++nc;
         na += cellPtr->nAtom();
         cellPtr = cellPtr->nextCellPtr();
      }
      TEST_ASSERT(na == locals.size());

      // Transform coordinates to Cartesian
      for (i = 0; i < nAtom; ++i) {
         for (j = 0; j < Dimension; ++j) {
            atoms[i].position()[j] *= lengths[j];
         }
      }
      cellList.update();

      // Find all neighbor pairs within a cutoff (use cell list)
      Cell::NeighborArray neighbors;
      CellAtom* cellAtomPtr1;
      CellAtom* cellAtomPtr2;
      Vector dr;
      int    nn;      // number of neighbors in a cell
      int    np = 0;  // Number of pairs within cutoff
      cellPtr = cellList.begin();
      while (cellPtr) {
         cellPtr->getNeighbors(neighbors);
         na = cellPtr->nAtom();
         ic = cellPtr->id();
         nn = neighbors.size();
         for (i = 0; i < na; ++i) {
            cellAtomPtr1 = neighbors[i];
            for (j = 0; j < na; ++j) {
               cellAtomPtr2 = neighbors[j];
               if (cellAtomPtr2 > cellAtomPtr1) {
                  dr.subtract(cellAtomPtr2->position(), cellAtomPtr1->position()); 
                  if (dr.square() < pairCutoffSq) {
                     ++np;
                  }
               }
            }
            for (j = na; j < nn; ++j) {
               cellAtomPtr2 = neighbors[j];
               dr.subtract(cellAtomPtr2->position(), cellAtomPtr1->position()); 
               if (dr.square() < pairCutoffSq) {
                  ++np;
               }
            }
         }
         cellPtr = cellPtr->nextCellPtr();
      }

      // Count neighbor pairs directly (N^2 loop)
      Atom* atomPtr1;
      Atom* atomPtr2;
      int nq = 0;
      for (i = 0; i < locals.size(); ++i) {
         atomPtr1 = &locals[i];
         for (j = 0; j < locals.size(); ++j) {
            atomPtr2 = &locals[j];
            if (atomPtr2 > atomPtr1) {
               dr.subtract(atomPtr2->position(), atomPtr1->position()); 
               if (dr.square() < pairCutoffSq) {
                  ++nq;
               }
            }
         }
         for (j = 0; j < ghosts.size(); ++j) {
            atomPtr2 = &ghosts[j];
            dr.subtract(atomPtr2->position(), atomPtr1->position()); 
            if (dr.square() < pairCutoffSq) {
               ++nq;
            }
         }
      }

      // Assert that number of neighobr pairs is same by either method.
      TEST_ASSERT(np == nq);
   }

   void testGetNeighbors2()
   {
      printMethod(TEST_FUNC);

      // Define cutoff
      double cutoff = 1.2;
      double pairCutoffSq = cutoff - 1.0E-10;
      pairCutoffSq = pairCutoffSq*pairCutoffSq;

      // Setup domain
      Vector lengths(4.0, 6.0, 8.0);
      Vector lower(2.0, 0.0, 4.0);
      Vector upper(4.0, 3.0, 8.0);
      Vector cutoffs;
      for (int i=0; i < 3; ++i) {
         lower[i] = lower[i]/lengths[i]; 
         upper[i] = upper[i]/lengths[i]; 
         cutoffs[i] = cutoff/lengths[i]; 
      }

      // Setup cell list grid
      const int nAtom = 200;
      int nCellCut = 3;
      cellList.allocate(nAtom, lower, upper, cutoffs, nCellCut);
      cellList.makeGrid(lower, upper, cutoffs, nCellCut);

      TEST_ASSERT(cellList.grid().dimension(0) == 11);
      TEST_ASSERT(cellList.grid().dimension(1) == 13);
      TEST_ASSERT(cellList.grid().dimension(2) == 16);
      TEST_ASSERT(eq(cellList.cellLength(0), 0.5/5.0));
      TEST_ASSERT(eq(cellList.cellLength(1), 0.5/7.0));
      TEST_ASSERT(eq(cellList.cellLength(2), 0.5/10.0));
      TEST_ASSERT(cellList.grid().size() == 11*13*16);

      Vector lowerGhost;
      Vector upperGhost;
      int i, j;
      for (i = 0; i < Dimension; ++i) {
         lowerGhost[i] = lower[i] - cutoffs[i];
         upperGhost[i] = upper[i] + cutoffs[i];
      }

      // Allocate Atom 
      AtomArray atoms;
      DPArray<Atom> locals;
      DPArray<Atom> ghosts;
      atoms.allocate(nAtom);
      locals.allocate(nAtom);
      ghosts.allocate(nAtom);

      // Create new random number generator
      Random random;
      random.setSeed(1298741);

      // Place atoms at random in extended region
      Vector pos;
      int ic;
      int nLocal = 0;
      int nGhost = 0;
      bool ghost;
      cellList.clear();
      for (i = 0; i < nAtom; ++i) {
         atoms[i].setTypeId(1);
         ghost = false;
         for (j = 0; j < Dimension; ++j) {
            pos[j] = random.uniform(lowerGhost[j], upperGhost[j]);
            if (pos[j] < lower[j]) 
               ghost = true;
            if (pos[j] > upper[j]) 
               ghost = true;
            TEST_ASSERT(pos[j] >= lowerGhost[j]);
            TEST_ASSERT(pos[j] <= upperGhost[j]);
         }
         atoms[i].position() = pos;
         if (ghost) {
            ++nGhost;
            ghosts.append(atoms[i]);
         } else {
            ++nLocal;
            locals.append(atoms[i]);
         }
         ic = cellList.cellIndexFromPosition(pos);
         TEST_ASSERT(ghost == cellList.cell(ic).isGhostCell());
         TEST_ASSERT(nGhost == ghosts.size());
         TEST_ASSERT(nLocal == locals.size());

         try {
            cellList.placeAtom(atoms[i]);
         } 
         catch (Exception e) {
            if (isIoProcessor()) {
               e.write(std::cout);
            }
            TEST_ASSERT(0);
         }

      }
      cellList.build();

      try { 
         cellList.isValid(); 
      } catch (Exception e) {
         TEST_ASSERT(0);
      }

      // Check that # atoms in local cells = # of local atoms
      int na = 0;
      int nc = 0;
      const Cell* cellPtr = cellList.begin();
      while (cellPtr) {
         ++nc;
         na += cellPtr->nAtom();
         cellPtr = cellPtr->nextCellPtr();
      }
      TEST_ASSERT(na == locals.size());

      // Transform coordinates to Cartesian
      for (i = 0; i < nAtom; ++i) {
         for (j = 0; j < Dimension; ++j) {
            atoms[i].position()[j] *= lengths[j];
         }
      }
      cellList.update();

      // Find all neighbor pairs within a cutoff (use cell list)
      Cell::NeighborArray neighbors;
      CellAtom* cellAtomPtr1;
      CellAtom* cellAtomPtr2;
      Vector dr;
      int    nn;      // number of neighbors in a cell
      int    np = 0;  // Number of pairs within cutoff
      cellPtr = cellList.begin();
      while (cellPtr) {
         cellPtr->getNeighbors(neighbors);
         na = cellPtr->nAtom();
         ic = cellPtr->id();
         nn = neighbors.size();
         for (i = 0; i < na; ++i) {
            cellAtomPtr1 = neighbors[i];
            for (j = 0; j < na; ++j) {
               cellAtomPtr2 = neighbors[j];
               if (cellAtomPtr2 > cellAtomPtr1) {
                  dr.subtract(cellAtomPtr2->position(), cellAtomPtr1->position()); 
                  if (dr.square() < pairCutoffSq) {
                     ++np;
                  }
               }
            }
            for (j = na; j < nn; ++j) {
               cellAtomPtr2 = neighbors[j];
               dr.subtract(cellAtomPtr2->position(), cellAtomPtr1->position()); 
               if (dr.square() < pairCutoffSq) {
                  ++np;
               }
            }
         }
         cellPtr = cellPtr->nextCellPtr();
      }

      // Count neighbor pairs directly (N^2 loop)
      Atom* atomPtr1;
      Atom* atomPtr2;
      int nq = 0;
      for (i = 0; i < locals.size(); ++i) {
         atomPtr1 = &locals[i];
         for (j = 0; j < locals.size(); ++j) {
            atomPtr2 = &locals[j];
            if (atomPtr2 > atomPtr1) {
               dr.subtract(atomPtr2->position(), atomPtr1->position()); 
               if (dr.square() < pairCutoffSq) {
                  ++nq;
               }
            }
         }
         for (j = 0; j < ghosts.size(); ++j) {
            atomPtr2 = &ghosts[j];
            dr.subtract(atomPtr2->position(), atomPtr1->position()); 
            if (dr.square() < pairCutoffSq) {
               ++nq;
            }
         }
      }

      // Assert that number of neighobr pairs is same by either method.
      TEST_ASSERT(np == nq);
   }
};

TEST_BEGIN(CellListTest)
TEST_ADD(CellListTest, testAllocate)
TEST_ADD(CellListTest, testAllocate2)
TEST_ADD(CellListTest, testMakeGrid)
TEST_ADD(CellListTest, testMakeGrid2)
TEST_ADD(CellListTest, testPlaceAtom)
TEST_ADD(CellListTest, testPlaceAtoms)
TEST_ADD(CellListTest, testAddRandomAtoms)
TEST_ADD(CellListTest, testAddRandomAtoms2)
TEST_ADD(CellListTest, testGetNeighbors)
TEST_ADD(CellListTest, testGetNeighbors2)
TEST_END(CellListTest)

#endif //ifndef DDMD_CELL_LIST_TEST_H
