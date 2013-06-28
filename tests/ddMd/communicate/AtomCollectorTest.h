#ifndef DDMD_ATOM_COLLECTOR_TEST_H
#define DDMD_ATOM_COLLECTOR_TEST_H

#include <ddMd/communicate/AtomCollector.h>
#include <ddMd/communicate/AtomDistributor.h>
#include <ddMd/communicate/Domain.h>
#include <ddMd/communicate/Buffer.h>
#include <ddMd/storage/AtomStorage.h>
#include <ddMd/storage/AtomIterator.h>
#include <ddMd/chemistry/Atom.h>
#include <util/space/Grid.h>
#include <util/mpi/MpiLogger.h>

#ifdef UTIL_MPI
#ifndef TEST_MPI
#define TEST_MPI
#endif
#endif

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>
#include <test/ParamFileTest.h>

using namespace Util;
using namespace DdMd;

class AtomCollectorTest: public ParamFileTest
{

   Boundary boundary;
   Domain domain;
   Buffer buffer;
   AtomStorage storage;
   AtomDistributor distributor;
   AtomCollector collector;
   int atomCount; // Number to be distributed by master
   int myRank;    // Communicator rank

public:

   virtual void setUp()
   {}

   void distribute()
   {
      std::ifstream atomposfile;

      // Set connections between objects
      domain.setBoundary(boundary);
      distributor.associate(domain, boundary, storage, buffer);

      #ifdef UTIL_MPI
      // Set communicators
      domain.setGridCommunicator(communicator());
      domain.setIoCommunicator(communicator());
      storage.setIoCommunicator(communicator());
      buffer.setIoCommunicator(communicator());
      distributor.setIoCommunicator(communicator());
      #else
      domain.setRank(0);
      #endif

      // Open parameter file
      #ifdef UTIL_MPI
      openFile("in/AtomDistributor.213");
      #else
      openFile("in/AtomDistributor.111");
      #endif

      // Reading parameter file and close
      domain.readParam(file());
      storage.readParam(file());
      buffer.readParam(file());
      distributor.readParam(file());
      closeFile();

      Vector boundarylength(6.0, 3.0, 9.0);
      boundary.setOrthorhombic(boundarylength);

      atomCount = 0; // Number to be distributed by master
      myRank = domain.gridRank();

      // If I am the master processor.
      if (myRank == 0) {
         int i;
         Vector  pos;
         Atom*   ptr;

         //atomposfile.open("in/Atompositions");
         openInputFile("in/Atompositions",atomposfile);
         atomposfile >> atomCount;

         #if UTIL_MPI
         // Initialize the sendbuffer.
         distributor.setup();
         #endif
         for(i = 0; i < atomCount; ++i) {
            ptr = distributor.newAtomPtr();
            ptr->setId(i);
            ptr->setTypeId(0);

            // Read a position from file.
            { 
               Vector r;
               atomposfile >> r;
               boundary.transformCartToGen(r, ptr->position());
            }

            //Use position vector for velocity
            ptr->velocity() = ptr->position();
            distributor.addAtom();
         }
         file().close();
         distributor.send();
      } else { // If I am not the master processor
         distributor.receive();
      }
   }

   void testDistribute()
   {
      printMethod(TEST_FUNC);

      distribute();

      int recvCount = storage.nAtom();
      AtomIterator iter;
      storage.begin(iter);
      for ( ; iter.notEnd(); ++iter) {
         TEST_ASSERT(domain.isInDomain(iter->position()));
      }

      #if 0
      #ifdef UTIL_MPI
      MpiLogger logger;
      logger.begin();
      std::cout << "Processor: " << myRank
                << ", recvCount = " << recvCount << std::endl;
      logger.end();
      #endif 
      #endif

      // Check that all atoms are accounted for after distribution.
      #ifdef UTIL_MPI
      int nRecvAll;
      communicator().Reduce(&recvCount, &nRecvAll, 1, MPI::INT, MPI::SUM, 0);
      if (myRank == 0) {
         //std::cout << "Total atom count = " << nRecvAll << std::endl;
         TEST_ASSERT(nRecvAll == atomCount);
      }
      #else
      //std::cout << "Total atom count = " << recvCount << std::endl;
      TEST_ASSERT(recvCount == atomCount);
      #endif

   }

   void testCollect1()
   {
      printMethod(TEST_FUNC);

      // Distribute atoms among processors
      distribute();

      // Collect atoms
      collector.associate(domain, storage, buffer);
      if (domain.isMaster()) {  
         collector.allocate(buffer.atomCapacity());
         collector.setup();
         Atom* atomPtr = collector.nextPtr();
         int i = 0;
         while (atomPtr) {
            //std::cout << atomPtr->id() 
            //          << "  " << atomPtr->position() << std::endl;
            atomPtr = collector.nextPtr();
            ++i;
         }
         TEST_ASSERT(i == atomCount);
      } else { 
         collector.send();
      }
   }

   void testCollect2()
   {
      printMethod(TEST_FUNC);

      // Distribute atoms among processors
      distribute();

      // Collect atoms
      collector.associate(domain, storage, buffer);
      if (domain.isMaster()) {  
         collector.allocate(6);
         collector.setup();
         Atom* atomPtr = collector.nextPtr();
         int i = 0;
         while (atomPtr) {
            //std::cout << atomPtr->id() 
            //          << "  " << atomPtr->position() << std::endl;
            atomPtr = collector.nextPtr();
            ++i;
         }
         TEST_ASSERT(i == atomCount);
      } else { 
         collector.send();
      }
   }

};

TEST_BEGIN(AtomCollectorTest)
TEST_ADD(AtomCollectorTest, testDistribute)
TEST_ADD(AtomCollectorTest, testCollect1)
TEST_ADD(AtomCollectorTest, testCollect2)
TEST_END(AtomCollectorTest)
#endif 
