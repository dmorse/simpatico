#ifndef DDMD_BOND_COLLECTOR_TEST_H
#define DDMD_BOND_COLLECTOR_TEST_H

#include <ddMd/communicate/AtomDistributor.h>
#include <ddMd/communicate/GroupDistributor.tpp>
#include <ddMd/communicate/GroupCollector.tpp>
#include <ddMd/communicate/Domain.h>
#include <ddMd/communicate/Buffer.h>
#include <ddMd/storage/AtomStorage.h>
#include <ddMd/storage/AtomIterator.h>
#include <ddMd/chemistry/Atom.h>
#include <ddMd/storage/BondStorage.h>
#include <util/param/Label.h>
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

class BondCollectorTest: public ParamFileTest
{

    Boundary boundary;
    Domain domain;
    Buffer buffer;
    AtomStorage atomStorage;
    BondStorage bondStorage;
    AtomDistributor atomDistributor;
    GroupDistributor<2> bondDistributor;
    GroupCollector<2> bondCollector;
    int atomCount; // Number of atoms to be distributed by master
    int bondCount; // Number to bonds be distributed by master

public:

   virtual void setUp()
   {}

   void distribute()
   {
      // Set connections between atomDistributors
      domain.setBoundary(boundary);
      atomDistributor.associate(domain, boundary, atomStorage, buffer);
      bondDistributor.associate(domain, atomStorage, bondStorage, buffer);

      #ifdef UTIL_MPI
      // Set communicators
      domain.setGridCommunicator(communicator());
      domain.setIoCommunicator(communicator());
      atomStorage.setIoCommunicator(communicator());
      bondStorage.setIoCommunicator(communicator());
      buffer.setIoCommunicator(communicator());
      atomDistributor.setIoCommunicator(communicator());
      bondDistributor.setIoCommunicator(communicator());
      #else
      domain.setRank(0);
      #endif

      // Open parameter file
      #ifdef UTIL_MPI
      openFile("in/GroupDistributor.213");
      #else
      openFile("in/GroupDistributor.111");
      #endif

      domain.readParam(file());
      buffer.readParam(file());
      atomStorage.readParam(file());
      atomDistributor.readParam(file());
      bondStorage.readParam(file());
      bondDistributor.readParam(file());

      // Finish reading parameter file
      closeFile();

      atomCount = 0; // Number of atoms to be distributed by master
      bondCount = 0; // Number to bonds be distributed by master
      int myRank = domain.gridRank();

      // If I am the master processor.
      std::ifstream configFile;
      if (myRank == 0) {
         //configFile.open("in/config");
         openInputFile("in/config",configFile);
         configFile >> Label("BOUNDARY");
         configFile >> boundary;
      }
      bcast(domain.communicator(), boundary, 0);

      if (myRank == 0) {

         // Read Max number of atoms to be distributed by the master processor
         configFile >> Label("ATOMS");
         configFile >> Label("nAtom") >> atomCount;

         //std::cout << std::endl;
         //std::cout << "Num Atoms to be distributed = " 
         //          << atomCount << std::endl;

         // Initialize the sendbuffer.
         atomDistributor.setup();

         // Fill the atom atomDistributors
         Atom*   ptr;
         int     id, typeId;
         for(int i = 0; i < atomCount; ++i) {
            ptr = atomDistributor.newAtomPtr();

            configFile >> id >> typeId;
            ptr->setId(id);
            ptr->setTypeId(typeId);

            // Read a position from file.
            if (UTIL_ORTHOGONAL) {
               configFile >> ptr->position();
            } else {
               Vector r;
               configFile >> r;
               boundary.transformCartToGen(r, ptr->position());
            }

            configFile >> ptr->velocity();
            atomDistributor.addAtom();

         }
         file().close();

         // Send any atoms not sent previously.
         atomDistributor.send();

      } else { // If I am not the master processor

         atomDistributor.receive();

      }

      int recvCount = atomStorage.nAtom();
      AtomIterator iter;
      atomStorage.begin(iter);
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

      // Read bonds
      if (myRank == 0) {

         // Read Max number of atoms to be distributed by the master processor
         configFile >> Label("BONDS");
         configFile >> Label("nBond") >> bondCount;

         //std::cout << std::endl;
         //std::cout << "Num Bonds to be distributed = " 
         //          << bondCount << std::endl;

         #if UTIL_MPI
         // Initialize the sendbuffer.
         bondDistributor.setup();
         #endif

         // Read bonds and add to bondDistributor
         int i;
         Group<2>* ptr;
         for(i = 0; i < bondCount; ++i) {
            ptr = bondDistributor.newPtr();
            configFile >> *ptr;
            bondDistributor.add();
         }
         TEST_ASSERT(i == bondCount);
         // std::cout << "# bonds distributed = " << i << std::endl;

         // Send any bonds not sent previously.
         bondDistributor.send();

      } else { // If I am not the master processor

         bondDistributor.receive();

      }

   }

   void testDistribute()
   {
      printMethod(TEST_FUNC);
      distribute();
   }

   void testCollect()
   {
      printMethod(TEST_FUNC);

      // Distribute groups among processors
      distribute();

      // Collect groups
      bondCollector.associate(domain, bondStorage, buffer);
      if (domain.isMaster()) {  
         bondCollector.allocate(20);
         bondCollector.setup();
         Group<2>* groupPtr = bondCollector.nextPtr();
         int i = 0;
         //std::cout << std::endl;
         while (groupPtr) {
            //std::cout << *groupPtr << std::endl;
            groupPtr = bondCollector.nextPtr();
            ++i;
         }
         //std::cout << "# bonds collected = " << i << std::endl;
         TEST_ASSERT(i == bondCount);
      } else { 
         bondCollector.send();
      }
   }

};

TEST_BEGIN(BondCollectorTest)
TEST_ADD(BondCollectorTest, testDistribute)
TEST_ADD(BondCollectorTest, testCollect)
TEST_END(BondCollectorTest)

#endif /* BOND_COLLECTOR_TEST_H */
