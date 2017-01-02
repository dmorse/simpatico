#ifndef MCMD_PAIR_LIST_TEST_H
#define MCMD_PAIR_LIST_TEST_H

#include <test/UnitTest.h>
#include <mcMd/neighbor/CellList.h>
#include <mcMd/neighbor/PairList.h>
#include <mcMd/neighbor/PairIterator.h>
#include <mcMd/chemistry/Atom.h>
#include <util/random/Random.h>

#include <iostream>
#include <fstream>

using namespace Util;
using namespace McMd;

class PairListTest : public UnitTest 
{

private:

   Boundary  boundary;
   PairList  pairList;

public:

   void setUp()
   {}

   void tearDown()
   {}

   void testInitialize()
   {
      printMethod(TEST_FUNC);
      double potentialCutoff = 1.2;
      int    nAtom           = 20;

      // Set up PairList
      Vector Lin(2.0, 3.0, 4.0);
      boundary.setOrthorhombic(Lin);  

      std::ifstream in;
      openInputFile("in/PairList", in);
      pairList.readParam(in);
      pairList.initialize(nAtom, potentialCutoff);

      try {
         pairList.isValid();
      } catch (Exception e) {
         std::cout << e.message();
         TEST_ASSERT(false);
      }

   }


   void testBuild()
   {
      printMethod(TEST_FUNC);
      const int    nAtom = 20;
      const double skin             = 0.2;
      const double potentialCutoff  = 1.2;

      Atom  *atom1Ptr, *atom2Ptr;
      double cutoff, rSq;
      int    i, j, nPair;

      cutoff = potentialCutoff + skin;

      // Initialize Boundary
      Vector Lin(2.0, 3.0, 4.0);
      boundary.setOrthorhombic(Lin);  

      // Initialize PairList
      std::ifstream in;
      openInputFile("in/PairList", in);
      pairList.readParam(in);
      pairList.initialize(nAtom, potentialCutoff);

      // Allocate Atoms and initialize Ids
      RArray<Atom>  atoms;
      Atom::allocate(nAtom, atoms);

      // Place Atoms at random
      Vector pos;
      Random random;
      random.setSeed(1098640);
      for (i=0; i < nAtom; ++i) {
         boundary.randomPosition(random, pos);
         atoms[i].setTypeId(1);
         atoms[i].position() = pos;
      }

      // Setup pair list (set up internal cell list)

      // Add all atoms to cell list
      pairList.setup(boundary);
      for (i=0; i < nAtom; ++i) {
         pairList.addAtom(atoms[i]);
      }

      // Build pairList from completed cell list
      pairList.build(boundary);

      // Use PairList::isValid() as test
      try {
         pairList.isValid();
      } catch (Exception e) {
         std::cout << e.message();
         TEST_ASSERT(0);
      }

      // Logical checks
      TEST_ASSERT(pairList.atomCapacity_ - 1 - pairList.tList1_ 
                     + pairList.nAtom1_ == pairList.nAtom_);
//    TEST_ASSERT(pairList.nAtom1_ - 1 == pairList.tList1_);
      TEST_ASSERT(pairList.first_[0] == 0);
      TEST_ASSERT(pairList.first_[pairList.nAtom1_] == pairList.nAtom2_);

      // Check that all pairs have distance less than cutoff^2
      for (PairIterator iter(pairList); iter.notEnd(); ++iter) {
         iter.getPair(atom1Ptr, atom2Ptr);
         rSq = boundary.distanceSq(atom1Ptr->position(), atom2Ptr->position());
         TEST_ASSERT(rSq < cutoff*cutoff);
      }


      // Check the total number of pairs
      nPair = 0;
      for (i = 0; i < nAtom; ++i) {
         for (j = i; j < nAtom; ++j) {
            rSq = boundary.distanceSq(atoms[i].position(), atoms[j].position());
            if (j > i) {
               if (rSq < cutoff*cutoff) {
                  ++nPair;
               }
            }
         }
      }
      TEST_ASSERT(nPair == pairList.nAtom2_);

      // Verbose output - list atom Ids for all pairs
      if (verbose() > 1) {
         for (PairIterator iter(pairList); iter.notEnd(); ++iter) {
            iter.getPair(atom1Ptr, atom2Ptr);
            printf(" %5i %5i \n", atom1Ptr->id(), atom2Ptr->id());
         }
         printf("nAtom1_ = %10i\n", pairList.nAtom1_);
         printf("tList1_ = %10i\n", pairList.tList1_);
         for (i = pairList.atomCapacity_ - 1; i > pairList.tList1_; --i) {
            printf(" %10i %10i \n", i, pairList.atom1Ptrs_[i]->id());
         }
      }

      Atom::deallocate();
   }

};


TEST_BEGIN(PairListTest)
TEST_ADD(PairListTest, testInitialize)
TEST_ADD(PairListTest, testBuild)
TEST_END(PairListTest)

#endif
