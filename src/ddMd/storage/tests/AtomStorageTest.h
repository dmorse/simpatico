#ifndef DDMD_ATOM_STORAGE_TEST_H
#define DDMD_ATOM_STORAGE_TEST_H

#ifdef UTIL_MPI
#define TEST_MPI
#endif

#include <ddMd/storage/AtomStorage.h>
#include <ddMd/storage/AtomIterator.h>
#include <ddMd/storage/GhostIterator.h>
#include <ddMd/chemistry/Atom.h>
#include <util/containers/DPArray.h>

#include <test/ParamFileTest.h>
#include <test/UnitTestRunner.h>

using namespace Util;
using namespace DdMd;

class AtomStorageTest : public ParamFileTest<AtomStorage>
{

public:

   virtual void setUp()
   { 
      #ifdef UTIL_MPI 
      object().setParamCommunicator(communicator());
      #endif

      openFile("in/AtomStorage"); 
      object().readParam(file()); 
   }

   void testReadParam();

   void testAddAtoms();

   void testAddRemoveAtoms();

   void testClearGhosts();

   void testIterators();

   void testSnapshot();

};

inline void AtomStorageTest::testReadParam()
{  
   printMethod(TEST_FUNC); 
   if (verbose() > 0) {
      object().writeParam(std::cout);
   }
}

inline void AtomStorageTest::testAddAtoms()
{
   printMethod(TEST_FUNC);

   Atom* ptr53 = object().addAtom(53);
   TEST_ASSERT(object().find(53) == ptr53);
   TEST_ASSERT(ptr53->id() == 53);
   TEST_ASSERT(object().nAtom() == 1);
   TEST_ASSERT(object().nGhost() == 0);
   TEST_ASSERT(object().isValid());

   Atom* ptr18 = object().addAtom(18);
   TEST_ASSERT(object().find(18) == ptr18);
   TEST_ASSERT(ptr18->id() == 18);
   TEST_ASSERT(object().find(53) == ptr53);
   TEST_ASSERT(ptr53->id() == 53);
   TEST_ASSERT(object().nAtom() == 2);
   TEST_ASSERT(object().nGhost() == 0);
   TEST_ASSERT(object().isValid());

   Atom* ptr35 = object().addGhost(35);
   TEST_ASSERT(object().find(53) == ptr53);
   TEST_ASSERT(ptr53->id() == 53);
   TEST_ASSERT(object().find(35) == ptr35);
   TEST_ASSERT(ptr35->id() == 35);
   TEST_ASSERT(object().find(18) == ptr18);
   TEST_ASSERT(ptr18->id() == 18);
   TEST_ASSERT(object().nAtom() == 2);
   TEST_ASSERT(object().nGhost() == 1);
   TEST_ASSERT(object().isValid());

}

void AtomStorageTest::testAddRemoveAtoms()
{
   printMethod(TEST_FUNC);

   // Add two atoms
   Atom* ptr53 = object().addAtom(53);
   Atom* ptr18 = object().addAtom(18);

   TEST_ASSERT(object().find(53) == ptr53);
   TEST_ASSERT(ptr53->id() == 53);
   TEST_ASSERT(object().find(18) == ptr18);
   TEST_ASSERT(ptr18->id() == 18);
   TEST_ASSERT(object().nAtom() == 2);
   TEST_ASSERT(object().nGhost() == 0);
   TEST_ASSERT(object().isValid());

   // Remove one atom
   object().removeAtom(ptr53);
   TEST_ASSERT(object().find(53) == 0);
   TEST_ASSERT(ptr53->id() < 0);
   TEST_ASSERT(object().find(18) == ptr18);
   TEST_ASSERT(ptr18->id() == 18);
   TEST_ASSERT(object().nAtom() == 1);
   TEST_ASSERT(object().nGhost() == 0);
   TEST_ASSERT(object().isValid());

   // Add two atoms
   Atom* ptr67 = object().addAtom(67);
   Atom* ptr44 = object().addAtom(44);
   TEST_ASSERT(object().find(53) == 0);
   TEST_ASSERT(object().find(18) == ptr18);
   TEST_ASSERT(ptr18->id() == 18);
   TEST_ASSERT(object().find(67) == ptr67);
   TEST_ASSERT(ptr67->id() == 67);
   TEST_ASSERT(object().find(44) == ptr44);
   TEST_ASSERT(ptr44->id() == 44);
   TEST_ASSERT(object().nAtom() == 3);
   TEST_ASSERT(object().nGhost() == 0);
   TEST_ASSERT(object().isValid());

   Atom* ptr35 = object().addGhost(35);
   TEST_ASSERT(object().find(35) == ptr35);
   TEST_ASSERT(ptr35->id() == 35);
   TEST_ASSERT(object().nGhost() == 1);
   Atom* ptr82 = object().addGhost(82);
   TEST_ASSERT(object().find(82) == ptr82);
   TEST_ASSERT(ptr82->id() == 82);
   TEST_ASSERT(object().find(35) == ptr35);
   TEST_ASSERT(ptr35->id() == 35);
   TEST_ASSERT(object().nAtom() == 3);
   TEST_ASSERT(object().nGhost() == 2);

   object().removeGhost(ptr35);
   TEST_ASSERT(object().find(53) == 0);
   TEST_ASSERT(object().find(35) == 0);
   TEST_ASSERT(ptr35->id() < 0);
   TEST_ASSERT(object().find(18) == ptr18);
   TEST_ASSERT(ptr18->id() == 18);
   TEST_ASSERT(object().find(67) == ptr67);
   TEST_ASSERT(ptr67->id() == 67);
   TEST_ASSERT(object().find(82) == ptr82);
   TEST_ASSERT(ptr82->id() == 82);
   TEST_ASSERT(object().find(44) == ptr44);
   TEST_ASSERT(ptr44->id() == 44);
   TEST_ASSERT(object().nAtom() == 3);
   TEST_ASSERT(object().nGhost() == 1);
   TEST_ASSERT(object().isValid());

}

void AtomStorageTest::testClearGhosts()
{
   printMethod(TEST_FUNC);

   DPArray<Atom> localAtoms;
   DPArray<Atom> ghostAtoms;

   localAtoms.allocate(object().atomCapacity());
   ghostAtoms.allocate(object().ghostCapacity());

   // Add atoms
   localAtoms.append(*object().addAtom(53));
   localAtoms.append(*object().addAtom(18));
   localAtoms.append(*object().addAtom(39));
   localAtoms.append(*object().addAtom(44));
   localAtoms.append(*object().addAtom(82));

   // Add ghosts
   ghostAtoms.append(*object().addGhost(35));
   ghostAtoms.append(*object().addGhost(84));
   ghostAtoms.append(*object().addGhost(17));
   ghostAtoms.append(*object().addGhost(94));

   TEST_ASSERT(object().nAtom() == 5);
   TEST_ASSERT(object().nGhost() == 4);
   TEST_ASSERT(object().isValid());

   object().clearGhosts();
   TEST_ASSERT(object().nAtom() == 5);
   TEST_ASSERT(object().nGhost() == 0);
   TEST_ASSERT(object().isValid());
}
 
void AtomStorageTest::testIterators()
{
   printMethod(TEST_FUNC);

   DPArray<Atom> localAtoms;
   DPArray<Atom> ghostAtoms;

   localAtoms.allocate(object().atomCapacity());
   ghostAtoms.allocate(object().ghostCapacity());

   // Add local atoms
   localAtoms.append(*object().addAtom(53));
   localAtoms.append(*object().addAtom(18));
   localAtoms.append(*object().addAtom(44));
   localAtoms.append(*object().addAtom(82));
   localAtoms.append(*object().addAtom(39));

   // Add ghosts
   ghostAtoms.append(*object().addGhost(35));
   ghostAtoms.append(*object().addGhost(17));
   TEST_ASSERT(object().isValid());
 
   AtomIterator localIter;
   int nLocal = 0; 
   for (object().begin(localIter); !localIter.atEnd(); ++localIter) {
      ++nLocal;
      //std::cout << localIter->id() << std::endl;
   }
   TEST_ASSERT(nLocal == object().nAtom());
   TEST_ASSERT(nLocal == 5);

   GhostIterator ghostIter;
   int nGhost = 0; 
   for (object().begin(ghostIter); !ghostIter.atEnd(); ++ghostIter) {
      ++nGhost;
      //std::cout << ghostIter->id() << std::endl;
   }
   TEST_ASSERT(nGhost == object().nGhost());
   TEST_ASSERT(nGhost == 2);

   object().removeAtom(&localAtoms[1]);
   --nLocal;
   TEST_ASSERT(object().isValid());
   TEST_ASSERT(nLocal == object().nAtom());
   TEST_ASSERT(nLocal == 4);

   object().removeGhost(&ghostAtoms[1]);
   --nGhost;
   TEST_ASSERT(object().isValid());
   TEST_ASSERT(nGhost == object().nGhost());
   TEST_ASSERT(nGhost == 1);

}

void AtomStorageTest::testSnapshot()
{
   printMethod(TEST_FUNC);

   DPArray<Atom> atoms;
   atoms.allocate(object().atomCapacity());

   // Add atoms
   atoms.append(*object().addAtom(53));
   atoms.append(*object().addAtom(35));
   atoms.append(*object().addAtom(18));
   atoms.append(*object().addAtom(44));
   atoms.append(*object().addAtom(17));
   TEST_ASSERT(object().isValid());

   Atom* ptr;

   ptr = object().find(53);
   TEST_ASSERT(ptr !=0 );
   ptr->position()[0] = 0.1;
   ptr->position()[1] = 1.1;
   ptr->position()[2] = 2.1;

   ptr = object().find(35);
   TEST_ASSERT(ptr !=0 );
   ptr->position()[0] = 0.2;
   ptr->position()[1] = 3.1;
   ptr->position()[2] = 2.8;

   ptr = object().find(18);
   TEST_ASSERT(ptr !=0 );
   ptr->position()[0] = 5.2;
   ptr->position()[1] = 2.1;
   ptr->position()[2] = 2.7;

   ptr = object().find(44);
   TEST_ASSERT(ptr !=0 );
   ptr->position()[0] = 5.2;
   ptr->position()[1] = 2.1;
   ptr->position()[2] = 2.7;

   ptr = object().find(17);
   TEST_ASSERT(ptr !=0 );
   ptr->position()[0] = 2.4;
   ptr->position()[1] = 8.8;
   ptr->position()[2] = 3.7;

   object().makeSnapshot();
   TEST_ASSERT( eq(object().maxSqDisplacement(), 0.0));
   object().clearSnapshot();

   object().makeSnapshot();
   ptr = object().find(44);
   ptr->position()[0] += -0.25;
   ptr->position()[1] +=  0.15;
   ptr = object().find(18);
   ptr->position()[0] +=  0.05;
   ptr->position()[1] +=  0.10;
   ptr->position()[2] +=  0.15;
   ptr = object().find(35);
   ptr->position()[0] +=  0.01;
   ptr->position()[1] +=  0.15;
   ptr->position()[2] +=  0.25;
   ptr = object().find(17); // maximum displacement
   ptr->position()[0] +=  0.20;
   ptr->position()[1] +=  0.05;
   ptr->position()[2] += -0.30;

   //std::cout << std::endl;
   //std::cout << object().maxSqDisplacement();

   TEST_ASSERT( eq(object().maxSqDisplacement(), 0.1325));
   object().clearSnapshot();

}
 
TEST_BEGIN(AtomStorageTest)
TEST_ADD(AtomStorageTest, testReadParam)
TEST_ADD(AtomStorageTest, testAddAtoms)
TEST_ADD(AtomStorageTest, testAddRemoveAtoms)
TEST_ADD(AtomStorageTest, testClearGhosts)
TEST_ADD(AtomStorageTest, testIterators)
TEST_ADD(AtomStorageTest, testSnapshot)
TEST_END(AtomStorageTest)

#endif
