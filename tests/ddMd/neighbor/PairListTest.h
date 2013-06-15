#ifndef DDMD_PAIR_LIST_TEST_H
#define DDMD_PAIR_LIST_TEST_H

#include <ddMd/neighbor/PairList.h>
#include <ddMd/neighbor/CellList.h>
#include <ddMd/neighbor/PairIterator.h>
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

class PairListTest : public UnitTest 
{

private:

   PairList      pairList;
   CellList      cellList;
   AtomArray     atoms;
   DPArray<Atom> locals;
   DPArray<Atom> ghosts;
   Vector        lower;
   Vector        upper;
   double        cutoff;
   double        cutoffSq;
   int           nAtom;
   int           pairCapacity;

public:

   PairListTest()
    : lower(0.0),
      upper(2.0, 3.0, 4.0),
      cutoff(1.2),
      cutoffSq(1.44),
      nAtom(200),
      pairCapacity(500)
   {}
   
   void setUp()
   {

      // Setup CellList
      cellList.allocate(nAtom, lower, upper, cutoff);
      cellList.makeGrid(lower, upper, cutoff);
      pairList.allocate(nAtom, pairCapacity, cutoff);

      #if 0
      TEST_ASSERT(cellList.grid().dimension(0) == 3);
      TEST_ASSERT(cellList.grid().dimension(1) == 4);
      TEST_ASSERT(cellList.grid().dimension(2) == 5);
      TEST_ASSERT(eq(cellList.cellLength(0), 2.0));
      TEST_ASSERT(eq(cellList.cellLength(1), 1.5));
      TEST_ASSERT(eq(cellList.cellLength(2), 4.0/3.0));
      TEST_ASSERT(cellList.grid().size() == 60);
      #endif

      // Allocate Atom 
      atoms.allocate(nAtom);
      locals.allocate(nAtom);
      ghosts.allocate(nAtom);

      for (int i = 0; i < nAtom; ++i) {
         atoms[i].setId(i);
      }

   }


   void tearDown()
   {}

   void makeConfiguration()
   {
      Vector lowerGhost;
      Vector upperGhost;
      int    i, j;
      for (i = 0; i < Dimension; ++i) {
          lowerGhost[i] = lower[i] - cutoff;
          upperGhost[i] = upper[i] + cutoff;
      }

      // Create new random number generator
      Random random;
      random.setSeed(1098640);

      // Place atoms at random in extended region
      Vector pos;
      int nLocal = 0;
      int nGhost = 0;
      int ic;
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
            e.write(std::cout);
            TEST_ASSERT(0);
         }

      }
      cellList.build();
      cellList.isValid(); 

      #if 0
      // Check that # atoms in local cells = # of local atoms
      int na = 0;
      const Cell* cellPtr = cellList.begin();
      while (cellPtr) {
         na += cellPtr->nAtom();
         cellPtr = cellPtr->nextCellPtr();
      }
      TEST_ASSERT(na == locals.size());
      #endif
   }

   void testCountNeighbors()
   {
      printMethod(TEST_FUNC);

      makeConfiguration();

      // Find all neighbors (cell list)
      Cell::NeighborArray neighbors;
      const Cell* cellPtr = cellList.begin();
      Atom* atom1Ptr;
      Atom* atom2Ptr;
      Vector dr;
      int    na;      // number of atoms in a cell
      int    nn;      // number of neighbors for a cell
      int    np = 0;  // Number of pairs within cutoff
      int    i, j;

      cellPtr = cellList.begin();
      while (cellPtr) {
         cellPtr->getNeighbors(neighbors);
         na = cellPtr->nAtom();
         nn = neighbors.size();
         for (i = 0; i < na; ++i) {
            atom1Ptr = neighbors[i];
            for (j = 0; j < na; ++j) {
               atom2Ptr = neighbors[j];
               if (atom2Ptr > atom1Ptr) {
                  dr.subtract(atom2Ptr->position(), atom1Ptr->position()); 
                  if (dr.square() <= cutoffSq) {
                     ++np;
                  }
               }
            }
            for (j = na; j < nn; ++j) {
               atom2Ptr = neighbors[j];
               dr.subtract(atom2Ptr->position(), atom1Ptr->position()); 
               if (dr.square() <= cutoffSq) {
                  ++np;
               }
            }
         }
         cellPtr = cellPtr->nextCellPtr();
      }
      //std::cout << "Total number of pairs = " << np << std::endl;

      // Count neighbor pairs (N^2 loop)
      int nq = 0;
      for (i = 0; i < locals.size(); ++i) {
         atom1Ptr = &locals[i];
         for (j = 0; j < locals.size(); ++j) {
            atom2Ptr = &locals[j];
            if (atom2Ptr > atom1Ptr) {
               dr.subtract(atom2Ptr->position(), atom1Ptr->position()); 
               if (dr.square() <= cutoffSq) {
                  ++nq;
               }
            }
         }
         for (j = 0; j < ghosts.size(); ++j) {
            atom2Ptr = &ghosts[j];
            dr.subtract(atom2Ptr->position(), atom1Ptr->position()); 
            if (dr.square() <= cutoffSq) {
               ++nq;
            }
         }
      }
      //std::cout << "Total number of pairs = " << nq << std::endl;

      TEST_ASSERT(np == nq);

      pairList.build(cellList);
      TEST_ASSERT(np == pairList.nPair());

   }

   void testPairIterator()
   {
      printMethod(TEST_FUNC);

      makeConfiguration();
      pairList.build(cellList);

      PairIterator iter;
      Vector       dr;
      double       drSq;
      Atom* atom1Ptr;
      Atom* atom2Ptr;
      int   i = 0;

      // Check that all pairs have distance less than cutoff^2
      pairList.begin(iter);
      for ( ; iter.notEnd(); ++iter) {
         iter.getPair(atom1Ptr, atom2Ptr);
         dr.subtract(atom1Ptr->position(), atom2Ptr->position());
         drSq = dr.square();
         
         TEST_ASSERT(drSq < cutoffSq);

         #if 0
         std::cout << Int(i,5 )
                   << Int(atom1Ptr->id(),5 ) << Int(atom2Ptr->id(),5 )
                   << std::endl;
         #endif

         ++i;
      }

   }

};

TEST_BEGIN(PairListTest)
TEST_ADD(PairListTest, testCountNeighbors)
TEST_ADD(PairListTest, testPairIterator)
TEST_END(PairListTest)

#endif
