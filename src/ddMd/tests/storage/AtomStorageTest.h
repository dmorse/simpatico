#ifndef DDMD_ATOM_STORAGE_TEST_H
#define DDMD_ATOM_STORAGE_TEST_H

#include <ddMd/storage/AtomStorage.h>
#include <ddMd/storage/AtomIterator.h>
#include <ddMd/storage/GhostIterator.h>
#include <ddMd/chemistry/Atom.h>
#include <util/containers/DPArray.h>
#include <util/containers/DArray.h>
#include <util/random/Random.h>

#ifdef UTIL_MPI
#ifndef TEST_MPI
#define TEST_MPI
#endif
#endif

#include <test/ParamFileTest.h>
#include <test/UnitTestRunner.h>

using namespace Util;
using namespace DdMd;

class AtomStorageTest : public ParamFileTest
{

private:

   AtomStorage storage_;

public:

   virtual void setUp()
   { 
      #ifdef UTIL_MPI 
      storage_.setIoCommunicator(communicator());
      #endif

      openFile("in/AtomStorage"); 
      storage_.readParam(file()); 
   }

   void testReadParam();

   void testAddAtoms();

   void testAddRemoveAtoms();

   void testClear();

   void testIterators();

   void testSnapshot();

   void testTransforms();

};

inline void AtomStorageTest::testReadParam()
{  
   printMethod(TEST_FUNC); 
   if (verbose() > 0) {
      storage_.writeParam(std::cout);
   }
}

inline void AtomStorageTest::testAddAtoms()
{
   printMethod(TEST_FUNC);

   TEST_ASSERT(storage_.atomReservoir_.size() == storage_.atomCapacity());
   TEST_ASSERT(storage_.atomSet_.size() == 0);
   TEST_ASSERT(storage_.nAtom() == 0);

   Atom* ptr53 = storage_.addAtom(53);
   TEST_ASSERT(storage_.find(53) == ptr53);
   TEST_ASSERT(ptr53->id() == 53);
   TEST_ASSERT(storage_.nAtom() == 1);
   TEST_ASSERT(storage_.nGhost() == 0);
   TEST_ASSERT(storage_.isValid());
   TEST_ASSERT(storage_.nAtom() == 1);
   TEST_ASSERT(storage_.atomReservoir_.size() == storage_.atomCapacity()-1);
   TEST_ASSERT(&storage_.atomSet_[0] == ptr53);
   TEST_ASSERT(ptr53 = &storage_.atoms_[0]);

   Atom* ptr18 = storage_.addAtom(18);
   TEST_ASSERT(storage_.find(18) == ptr18);
   TEST_ASSERT(ptr18->id() == 18);
   TEST_ASSERT(storage_.find(53) == ptr53);
   TEST_ASSERT(ptr53->id() == 53);
   TEST_ASSERT(storage_.nAtom() == 2);
   TEST_ASSERT(storage_.atomReservoir_.size() == storage_.atomCapacity()-2);
   TEST_ASSERT(&storage_.atomSet_[0] == ptr53);
   TEST_ASSERT(&storage_.atomSet_[1] == ptr18);
   TEST_ASSERT(ptr53 = &storage_.atoms_[0]);
   TEST_ASSERT(ptr18 = &storage_.atoms_[1]);
   TEST_ASSERT(storage_.nGhost() == 0);
   TEST_ASSERT(storage_.isValid());

   Atom* ptr35 = storage_.addGhost(35);
   TEST_ASSERT(storage_.find(53) == ptr53);
   TEST_ASSERT(ptr53->id() == 53);
   TEST_ASSERT(storage_.find(35) == ptr35);
   TEST_ASSERT(ptr35->id() == 35);
   TEST_ASSERT(storage_.find(18) == ptr18);
   TEST_ASSERT(ptr18->id() == 18);
   TEST_ASSERT(storage_.nAtom() == 2);
   TEST_ASSERT(storage_.nGhost() == 1);
   TEST_ASSERT(storage_.isValid());

}

void AtomStorageTest::testAddRemoveAtoms()
{
   printMethod(TEST_FUNC);

   // Add two atoms
   Atom* ptr53 = storage_.addAtom(53);
   Atom* ptr18 = storage_.addAtom(18);

   TEST_ASSERT(storage_.find(53) == ptr53);
   TEST_ASSERT(ptr53->id() == 53);
   TEST_ASSERT(storage_.find(18) == ptr18);
   TEST_ASSERT(ptr18->id() == 18);
   TEST_ASSERT(storage_.nAtom() == 2);
   TEST_ASSERT(storage_.nGhost() == 0);
   TEST_ASSERT(storage_.isValid());

   // Remove one atom
   storage_.removeAtom(ptr53);
   TEST_ASSERT(storage_.find(53) == 0);
   TEST_ASSERT(ptr53->id() < 0);
   TEST_ASSERT(storage_.find(18) == ptr18);
   TEST_ASSERT(ptr18->id() == 18);
   TEST_ASSERT(storage_.nAtom() == 1);
   TEST_ASSERT(storage_.nGhost() == 0);
   TEST_ASSERT(storage_.isValid());

   // Add two atoms
   Atom* ptr67 = storage_.addAtom(67);
   Atom* ptr44 = storage_.addAtom(44);
   TEST_ASSERT(storage_.find(53) == 0);
   TEST_ASSERT(storage_.find(18) == ptr18);
   TEST_ASSERT(ptr18->id() == 18);
   TEST_ASSERT(storage_.find(67) == ptr67);
   TEST_ASSERT(ptr67->id() == 67);
   TEST_ASSERT(storage_.find(44) == ptr44);
   TEST_ASSERT(ptr44->id() == 44);
   TEST_ASSERT(storage_.nAtom() == 3);
   TEST_ASSERT(storage_.nGhost() == 0);
   TEST_ASSERT(storage_.isValid());

   Atom* ptr35 = storage_.addGhost(35);
   TEST_ASSERT(storage_.find(35) == ptr35);
   TEST_ASSERT(ptr35->id() == 35);
   TEST_ASSERT(storage_.nGhost() == 1);
   Atom* ptr82 = storage_.addGhost(82);
   TEST_ASSERT(storage_.find(82) == ptr82);
   TEST_ASSERT(ptr82->id() == 82);
   TEST_ASSERT(storage_.find(35) == ptr35);
   TEST_ASSERT(ptr35->id() == 35);
   TEST_ASSERT(storage_.nAtom() == 3);
   TEST_ASSERT(storage_.nGhost() == 2);

   storage_.removeGhost(ptr35);
   TEST_ASSERT(storage_.find(53) == 0);
   TEST_ASSERT(storage_.find(35) == 0);
   TEST_ASSERT(ptr35->id() < 0);
   TEST_ASSERT(storage_.find(18) == ptr18);
   TEST_ASSERT(ptr18->id() == 18);
   TEST_ASSERT(storage_.find(67) == ptr67);
   TEST_ASSERT(ptr67->id() == 67);
   TEST_ASSERT(storage_.find(82) == ptr82);
   TEST_ASSERT(ptr82->id() == 82);
   TEST_ASSERT(storage_.find(44) == ptr44);
   TEST_ASSERT(ptr44->id() == 44);
   TEST_ASSERT(storage_.nAtom() == 3);
   TEST_ASSERT(storage_.nGhost() == 1);
   TEST_ASSERT(storage_.isValid());
}

void AtomStorageTest::testClear()
{
   printMethod(TEST_FUNC);

   //DPArray<Atom> localAtoms;
   DPArray<Atom> ghostAtoms;
   //localAtoms.allocate(storage_.atomCapacity());
   ghostAtoms.allocate(storage_.ghostCapacity());

   // Add atoms
   Atom* ptr53 = storage_.addAtom(53);
   TEST_ASSERT(storage_.find(53) == ptr53);
   TEST_ASSERT(ptr53->id() == 53);

   Atom* ptr18 = storage_.addAtom(18);
   TEST_ASSERT(ptr18 != ptr53);
   TEST_ASSERT(storage_.find(53) == ptr53);
   TEST_ASSERT(ptr53->id() == 53);
   TEST_ASSERT(storage_.find(18) == ptr18);
   TEST_ASSERT(ptr18->id() == 18);

   Atom* ptr39 = storage_.addAtom(39);
   TEST_ASSERT(ptr39 != ptr53);
   TEST_ASSERT(ptr39 != ptr18);
   TEST_ASSERT(storage_.find(53) == ptr53);
   TEST_ASSERT(ptr53->id() == 53);
   TEST_ASSERT(storage_.find(18) == ptr18);
   TEST_ASSERT(ptr18->id() == 18);
   TEST_ASSERT(storage_.find(39) == ptr39);
   TEST_ASSERT(ptr39->id() == 39);

   Atom* ptr44 = storage_.addAtom(44);
   TEST_ASSERT(ptr44 != ptr53);
   TEST_ASSERT(ptr44 != ptr18);
   TEST_ASSERT(ptr44 != ptr39);
   TEST_ASSERT(ptr39 != ptr53);
   TEST_ASSERT(ptr39 != ptr18);
   TEST_ASSERT(storage_.find(53) == ptr53);
   TEST_ASSERT(ptr53->id() == 53);

   Atom* ptr82 = storage_.addAtom(82);
   TEST_ASSERT(ptr82 != ptr53);
   TEST_ASSERT(ptr82 != ptr44);
   TEST_ASSERT(ptr82 != ptr18);
   TEST_ASSERT(ptr82 != ptr39);
   TEST_ASSERT(storage_.find(82) == ptr82);
   TEST_ASSERT(ptr82->id() == 82);
   TEST_ASSERT(storage_.find(53) == ptr53);
   TEST_ASSERT(ptr53->id() == 53);
   
   //localAtoms.append(*(storage_.addAtom(53)));
   //localAtoms.append(*(storage_.addAtom(18)));
   //localAtoms.append(*(storage_.addAtom(39)));
   //localAtoms.append(*(storage_.addAtom(44)));
   //localAtoms.append(*(storage_.addAtom(82)));
   
   TEST_ASSERT(storage_.isValid());
   TEST_ASSERT(storage_.nAtom() == 5);
   TEST_ASSERT(storage_.nGhost() == 0);

   // Add ghosts
   ghostAtoms.append(*storage_.addGhost(35));
   ghostAtoms.append(*storage_.addGhost(84));
   ghostAtoms.append(*storage_.addGhost(17));
   ghostAtoms.append(*storage_.addGhost(94));
   TEST_ASSERT(storage_.isValid());
   TEST_ASSERT(storage_.nAtom() == 5);
   TEST_ASSERT(storage_.nGhost() == 4);

   storage_.clearGhosts();
   TEST_ASSERT(storage_.isValid());
   TEST_ASSERT(storage_.nAtom() == 5);
   TEST_ASSERT(storage_.nGhost() == 0);

   storage_.clearAtoms();
   TEST_ASSERT(storage_.isValid());
   TEST_ASSERT(storage_.nAtom() == 0);
   TEST_ASSERT(storage_.nGhost() == 0);

}
 
void AtomStorageTest::testIterators()
{
   printMethod(TEST_FUNC);

   DPArray<Atom> localAtoms;
   DPArray<Atom> ghostAtoms;

   localAtoms.allocate(storage_.atomCapacity());
   ghostAtoms.allocate(storage_.ghostCapacity());

   // Add local atoms
   localAtoms.append(*storage_.addAtom(53));
   localAtoms.append(*storage_.addAtom(18));
   localAtoms.append(*storage_.addAtom(44));
   localAtoms.append(*storage_.addAtom(82));
   localAtoms.append(*storage_.addAtom(39));

   // Add ghosts
   ghostAtoms.append(*storage_.addGhost(35));
   ghostAtoms.append(*storage_.addGhost(17));
   TEST_ASSERT(storage_.isValid());
 
   AtomIterator localIter;
   int nLocal = 0; 
   for (storage_.begin(localIter); localIter.notEnd(); ++localIter) {
      ++nLocal;
      //std::cout << localIter->id() << std::endl;
   }
   TEST_ASSERT(nLocal == storage_.nAtom());
   TEST_ASSERT(nLocal == 5);

   GhostIterator ghostIter;
   int nGhost = 0; 
   for (storage_.begin(ghostIter); ghostIter.notEnd(); ++ghostIter) {
      ++nGhost;
      //std::cout << ghostIter->id() << std::endl;
   }
   TEST_ASSERT(nGhost == storage_.nGhost());
   TEST_ASSERT(nGhost == 2);

   storage_.removeAtom(&localAtoms[1]);
   --nLocal;
   TEST_ASSERT(storage_.isValid());
   TEST_ASSERT(nLocal == storage_.nAtom());
   TEST_ASSERT(nLocal == 4);

   storage_.removeGhost(&ghostAtoms[1]);
   --nGhost;
   TEST_ASSERT(storage_.isValid());
   TEST_ASSERT(nGhost == storage_.nGhost());
   TEST_ASSERT(nGhost == 1);

}

void AtomStorageTest::testSnapshot()
{
   printMethod(TEST_FUNC);

   DPArray<Atom> atoms;
   atoms.allocate(storage_.atomCapacity());

   // Add atoms
   atoms.append(*storage_.addAtom(53));
   atoms.append(*storage_.addAtom(35));
   atoms.append(*storage_.addAtom(18));
   atoms.append(*storage_.addAtom(44));
   atoms.append(*storage_.addAtom(17));
   TEST_ASSERT(storage_.isValid());

   Atom* ptr;

   ptr = storage_.find(53);
   TEST_ASSERT(ptr !=0 );
   ptr->position()[0] = 0.1;
   ptr->position()[1] = 1.1;
   ptr->position()[2] = 2.1;

   ptr = storage_.find(35);
   TEST_ASSERT(ptr !=0 );
   ptr->position()[0] = 0.2;
   ptr->position()[1] = 3.1;
   ptr->position()[2] = 2.8;

   ptr = storage_.find(18);
   TEST_ASSERT(ptr !=0 );
   ptr->position()[0] = 5.2;
   ptr->position()[1] = 2.1;
   ptr->position()[2] = 2.7;

   ptr = storage_.find(44);
   TEST_ASSERT(ptr !=0 );
   ptr->position()[0] = 5.2;
   ptr->position()[1] = 2.1;
   ptr->position()[2] = 2.7;

   ptr = storage_.find(17);
   TEST_ASSERT(ptr !=0 );
   ptr->position()[0] = 2.4;
   ptr->position()[1] = 8.8;
   ptr->position()[2] = 3.7;

   
   OrthorhombicBoundary boundary;
   Vector L(2.0, 5.0, 4.0);
   boundary.setOrthorhombic(L);

   // Transform too and from generalized coordinates.
   // This simply marks the coordinate system as Cartesian.
   AtomIterator iter;
   Vector Rg;
   for (storage_.begin(iter); iter.notEnd(); ++iter) {
      boundary.transformCartToGen(iter->position(), Rg);
      iter->position() = Rg; 
   }
   storage_.transformGenToCart(boundary);

   storage_.makeSnapshot();
   TEST_ASSERT( eq(storage_.maxSqDisplacement(), 0.0));
   storage_.clearSnapshot();

   storage_.makeSnapshot();
   ptr = storage_.find(44);
   ptr->position()[0] += -0.25;
   ptr->position()[1] +=  0.15;
   ptr = storage_.find(18);
   ptr->position()[0] +=  0.05;
   ptr->position()[1] +=  0.10;
   ptr->position()[2] +=  0.15;
   ptr = storage_.find(35);
   ptr->position()[0] +=  0.01;
   ptr->position()[1] +=  0.15;
   ptr->position()[2] +=  0.25;
   ptr = storage_.find(17); // maximum displacement
   ptr->position()[0] +=  0.20;
   ptr->position()[1] +=  0.05;
   ptr->position()[2] += -0.30;

   //std::cout << std::endl;
   //std::cout << storage_.maxSqDisplacement();

   TEST_ASSERT(eq(storage_.maxSqDisplacement(), 0.1325));
   storage_.clearSnapshot();

}
 
void AtomStorageTest::testTransforms()
{
   printMethod(TEST_FUNC);

   DPArray<Atom>  localAtoms;
   DPArray<Atom>  ghostAtoms;
   DArray<Vector> gPositions;
   DArray<Vector> cPositions;
   Random        random;

   localAtoms.allocate(storage_.atomCapacity());
   ghostAtoms.allocate(storage_.ghostCapacity());
   int totalCapacity = storage_.atomCapacity() + storage_.ghostCapacity();
   gPositions.allocate(totalCapacity); 
   cPositions.allocate(totalCapacity);
   random.setSeed(274454136);

   // Add local atoms
   localAtoms.append(*storage_.addAtom(53));
   localAtoms.append(*storage_.addAtom(18));
   localAtoms.append(*storage_.addAtom(44));
   localAtoms.append(*storage_.addAtom(82));
   localAtoms.append(*storage_.addAtom(39));
   localAtoms.append(*storage_.addAtom(21));
   localAtoms.append(*storage_.addAtom(76));

   // Add ghosts
   ghostAtoms.append(*storage_.addGhost(35));
   ghostAtoms.append(*storage_.addGhost(17));
   ghostAtoms.append(*storage_.addGhost(92));
   ghostAtoms.append(*storage_.addGhost(28));
   ghostAtoms.append(*storage_.addGhost(73));

   TEST_ASSERT(storage_.isValid());
 
   AtomIterator localIter;
   int nLocal = 0; 
   int i = 0;
   int j;
   for (storage_.begin(localIter); localIter.notEnd(); ++localIter) {
      for (j = 0; j < Dimension; ++j) {
         localIter->position()[j] = random.uniform(-0.5, 1.5);
      }
      gPositions[i] = localIter->position();
      ++nLocal;
      ++i;
   }
   TEST_ASSERT(nLocal == storage_.nAtom());
   TEST_ASSERT(nLocal == 7);

   GhostIterator ghostIter;
   int nGhost = 0; 
   for (storage_.begin(ghostIter); ghostIter.notEnd(); ++ghostIter) {
      for (j = 0; j < Dimension; ++j) {
         ghostIter->position()[j] = random.uniform(-0.5, 1.5);
      }
      gPositions[i] = ghostIter->position();
      ++nGhost;
      ++i;
   }
   TEST_ASSERT(nGhost == storage_.nGhost());
   TEST_ASSERT(nGhost == 5);

   OrthorhombicBoundary boundary;
   Vector L(2.0, 3.0, 4.0);
   boundary.setOrthorhombic(L);

   TEST_ASSERT(!storage_.isCartesian());
   storage_.transformGenToCart(boundary);
   TEST_ASSERT(storage_.isCartesian());
   storage_.transformCartToGen(boundary);
   TEST_ASSERT(!storage_.isCartesian());

   // Check against stored positions
   i = 0;
   for (storage_.begin(localIter); localIter.notEnd(); ++localIter) {
      for (j = 0; j < Dimension; ++j) {
         TEST_ASSERT(eq(localIter->position()[j], gPositions[i][j]));
      }
      ++i;
   }
   for (storage_.begin(ghostIter); ghostIter.notEnd(); ++ghostIter) {
      for (j = 0; j < Dimension; ++j) {
         TEST_ASSERT(eq(ghostIter->position()[j], gPositions[i][j]));
      }
      ++i;
   }

}

TEST_BEGIN(AtomStorageTest)
TEST_ADD(AtomStorageTest, testReadParam)
TEST_ADD(AtomStorageTest, testAddAtoms)
TEST_ADD(AtomStorageTest, testAddRemoveAtoms)
TEST_ADD(AtomStorageTest, testClear)
TEST_ADD(AtomStorageTest, testIterators)
TEST_ADD(AtomStorageTest, testSnapshot)
TEST_ADD(AtomStorageTest, testTransforms)
TEST_END(AtomStorageTest)

#endif
