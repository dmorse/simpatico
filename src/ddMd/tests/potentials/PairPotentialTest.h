#ifndef DDMD_PAIR_POTENTIAL_TEST_H
#define DDMD_PAIR_POTENTIAL_TEST_H

#include <ddMd/potentials/pair/PairPotential.h>
#include <ddMd/potentials/pair/PairPotentialImpl.h>
#include <ddMd/storage/AtomStorage.h>
#include <ddMd/storage/AtomIterator.h>
#include <ddMd/storage/GhostIterator.h>
#include <ddMd/communicate/Domain.h>
#include <ddMd/chemistry/Atom.h>
#include <simp/interaction/pair/DpdPair.h>
#include <util/boundary/Boundary.h>
#include <util/random/Random.h>

#ifdef UTIL_MPI
#ifndef TEST_MPI
#define TEST_MPI
#endif
#endif

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>
#include <test/ParamFileTest.h>

using namespace Util;
using namespace Simp;
using namespace DdMd;

class PairPotentialTest: public ParamFileTest
{

private:

      Boundary boundary;
      Domain   domain;
      AtomStorage storage;
      PairPotentialImpl<DpdPair> pairPotential;

public:

   virtual void setUp()
   {
      pairPotential.setNAtomType(1);
      pairPotential.associate(domain, boundary, storage);
      domain.setBoundary(boundary);

      #ifdef UTIL_MPI
      domain.setGridCommunicator(communicator());
      domain.setIoCommunicator(communicator());
      storage.setIoCommunicator(communicator());
      pairPotential.setIoCommunicator(communicator());
      #endif

      // Read parameter file
      openFile("in/PairPotential");
      domain.readParam(file());
      storage.readParam(file());
      pairPotential.readParam(file());
      closeFile();
   }

   void readAtoms(const char* filename)
   {

      int atomCount = 0; // Number to be distributed by master
      int i;
      Vector  pos;
      Atom*   ptr;

      std::ifstream atomposfile;
      atomposfile.open(filename);

      // Read Max # of atoms to be distributed.
      atomposfile >> atomCount;

      //std::cout << std::endl;
      //std::cout << "Num Atoms to be distributed = " << atomCount << std::endl;

      // Read atoms 
      for(i = 0; i < atomCount; ++i) {
         ptr = storage.newAtomPtr();
         ptr->setId(i);
         ptr->setTypeId(0);

         //Read a position from file.
         atomposfile >> ptr->position();

         //Use position vector for velocity for now
         ptr->velocity() = ptr->position();

         storage.addNewAtom();

      }
      file().close();

   }

   void randomAtoms(int nAtom, const Vector& lower, const Vector& upper, 
                    double cutoff)
   {

      // Set coordinate bounds
      Vector  lowerGhost;
      Vector  upperGhost;
      Vector  length;
      int    i, j;
      for (i = 0; i < Dimension; ++i) {
          length[i] = upper[i] - lower[i];
          lowerGhost[i] = (lower[i] - cutoff)/length[i];
          upperGhost[i] = (upper[i] + cutoff)/length[i];
      }

      // Create new random number generator
      Random random;
      random.setSeed(20);

      // Place atoms at random in extended region
      Vector pos;
      Atom*  ptr;
      bool   ghost;
      storage.clearAtoms();
      for (i = 0; i < nAtom; ++i) {
         ghost = false;
         for (j = 0; j < Dimension; ++j) {
            pos[j] = random.uniform(lowerGhost[j], upperGhost[j]);
            if (pos[j] < lower[j]) ghost = true;
            if (pos[j] > upper[j]) ghost = true;
            TEST_ASSERT(pos[j] >= lowerGhost[j]);
            TEST_ASSERT(pos[j] <= upperGhost[j]);
         }
         if (ghost) {
            ptr = storage.newGhostPtr();
            ptr->setId(i);
            ptr->setTypeId(0);
            ptr->position() = pos;
            storage.addNewGhost();
         } else {
            ptr = storage.newAtomPtr();
            ptr->setId(i);
            ptr->setTypeId(0);
            ptr->position() = pos;
            storage.addNewAtom();
         }
         std::cout << pos << "   " << ghost << std::endl;
      }
      TEST_ASSERT(!storage.isCartesian());
   }

   void zeroForces()
   {
      AtomIterator  atomIter;
      GhostIterator ghostIter;

      storage.begin(atomIter);
      for ( ; atomIter.notEnd(); ++atomIter) {
         atomIter->force().zero();
      }
      storage.begin(ghostIter);
      for ( ; ghostIter.notEnd(); ++ghostIter) {
         ghostIter->force().zero();
      }

   }

   void writeForces()
   {
      std::cout << std::endl;

      AtomIterator  atomIter;
      storage.begin(atomIter);
      for ( ; atomIter.notEnd(); ++atomIter) {
         std::cout << atomIter->force() << "  " << 0 << std::endl;
      }

      #if 0
      GhostIterator ghostIter;
      storage.begin(ghostIter);
      for ( ; ghostIter.notEnd(); ++ghostIter) {
         std::cout << ghostIter->force() << "  " << 1 << std::endl;
      }
      #endif

   }

   #if 0
   void saveForces()
   {
      AtomIterator  atomIter;
      GhostIterator ghostIter;

      storage.begin(atomIter);
      for ( ; atomIter.notEnd(); ++atomIter) {
         force_[i] = atomIter->force();
      }
      storage.begin(ghostIter);
      for ( ; ghostIter.notEnd(); ++ghostIter) {
         force_[i] = ghostIter->force();
      }

   }
   #endif

   void testRead1()
   {
      printMethod(TEST_FUNC);
      readAtoms("in/positions1");
     
      TEST_ASSERT(storage.nAtom()  == 3); 
      TEST_ASSERT(storage.nGhost() == 0); 

      AtomIterator iter;
      std::cout << std::endl;
      storage.begin(iter);
      for ( ; iter.notEnd(); ++iter) {
         std::cout << iter->position() << std::endl;
      }      

      pairPotential.setMethodId(2); // N^2 loop
      pairPotential.computeForces();

      std::cout << std::endl;
      storage.begin(iter);
      for ( ; iter.notEnd(); ++iter) {
         std::cout << iter->force() << std::endl;
      }      
   }

   void testRandom1()
   {
      printMethod(TEST_FUNC);

      const int nAtom = 120;
      double cutoff   = 1.2;
      Vector lower(0.0);
      Vector upper(2.0, 3.0, 4.0);

      boundary.setOrthorhombic(upper);
      randomAtoms(nAtom, lower, upper, cutoff);

      TEST_ASSERT(!storage.isCartesian());
      pairPotential.buildCellList();
      storage.transformGenToCart(boundary);
      pairPotential.buildPairList();

      zeroForces();
      pairPotential.setMethodId(0);
      pairPotential.computeForces();
      writeForces();
 
      zeroForces();
      pairPotential.setMethodId(1);
      pairPotential.computeForces();
      writeForces();

      zeroForces();
      pairPotential.setMethodId(2);
      pairPotential.computeForces();
      writeForces();

   }

};

TEST_BEGIN(PairPotentialTest)
TEST_ADD(PairPotentialTest, testRead1)
TEST_ADD(PairPotentialTest, testRandom1)
TEST_END(PairPotentialTest)

#endif 
