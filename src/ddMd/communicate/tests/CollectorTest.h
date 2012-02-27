#ifndef COLLECTOR_TEST_H
#define COLLECTOR_TEST_H

#include <ddMd/communicate/Collector.h>
#include <ddMd/communicate/AtomDistributor.h>
#include <ddMd/communicate/Domain.h>
#include <ddMd/communicate/Buffer.h>
#include <ddMd/storage/AtomStorage.h>
#include <ddMd/storage/AtomIterator.h>
#include <ddMd/chemistry/Atom.h>
#include <util/space/Grid.h>
#include <util/mpi/MpiLogger.h>

#ifdef  UTIL_MPI
   #ifndef TEST_MPI
      #define TEST_MPI
   #endif
#endif

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>
#include <test/ParamFileTest.h>

using namespace Util;
using namespace DdMd;

class CollectorTest: public ParamFileTest<Collector>
{
public:

   virtual void setUp()
   {}

   void testCollect()
   {
      printMethod(TEST_FUNC);

      Boundary         boundary;
      Domain           domain;
      Buffer           buffer;
      AtomStorage      storage;
      AtomDistributor  distributor;
      std::ifstream    atomposfile;

      // Set connections between objects
      domain.setBoundary(boundary);
      distributor.associate(domain, boundary, storage, buffer);

      #ifdef UTIL_MPI
      // Set communicators
      domain.setGridCommunicator(communicator());
      domain.setParamCommunicator(communicator());
      storage.setParamCommunicator(communicator());
      buffer.setParamCommunicator(communicator());
      distributor.setParamCommunicator(communicator());
      object().setParamCommunicator(communicator());
      #else
      domain.setRank(0);
      #endif

      // Open parameter file
      #ifdef UTIL_MPI
         openFile("in/Collector.213");
      #else
         openFile("in/Collector.111");
      #endif

      domain.readParam(file());

      // Allocate/initialize AtomStorage
      storage.readParam(file());

      // Initialize buffer
      buffer.readParam(file());

      // Initialize AtomDistributor object
      distributor.readParam(file());

      // Finish reading parameter file
      closeFile();

      Vector boundarylength(6.0, 3.0, 9.0);
      boundary.setLengths(boundarylength);

      int atomCount = 0; // Number to be distributed by master
      int myRank    = domain.gridRank();

      // If I am the master processor.
      if (myRank == 0) {
         int i;
         Vector  pos;
         Atom*   ptr;

         atomposfile.open("in/Atompositions");
         // Read Max number of atoms to be distributed by the master processor
         atomposfile >> atomCount;

         std::cout << std::endl;
         std::cout << "Num Atoms to be distributed = " << atomCount << std::endl;

         #if UTIL_MPI
         // Initialize the sendbuffer.
         distributor.initSendBuffer();
         #endif

         // Fill the atom objects
         for(i = 0; i < atomCount; ++i) {
            ptr = distributor.newAtomPtr();
            ptr->setId(i);
            ptr->setTypeId(0);

            //Assign a random position within the boundary for the atom.
            atomposfile >> ptr->position();

            //Use position vector for velocity for now
            ptr->velocity() = ptr->position();

            distributor.addAtom();

         }
         file().close();

         // Send any atoms not sent previously.
         distributor.send();

      } else { // If I am not the master processor

         distributor.receive();

      }

      int recvCount = storage.nAtom();
      AtomIterator iter;
      storage.begin(iter);
      for ( ; !iter.atEnd(); ++iter) {
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

       // If master processor
       if (domain.isMaster()) {  
   
          object().initialize(storage, domain, buffer);
          Atom* atomPtr = object().nextPtr();
          int i = 0;
          while (atomPtr) {
             std::cout << i << "  " << atomPtr->position() << std::endl;
             atomPtr = object().nextPtr();
             ++i;
          }
   
       } else { // if not master processor
   
          object().send(storage, domain, buffer);
   
       }
   }
};

TEST_BEGIN(CollectorTest)
TEST_ADD(CollectorTest, testCollect)
TEST_END(CollectorTest)

#endif /* COLLECTOR_TEST_H */
