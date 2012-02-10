#ifndef BOND_DISTRIBUTOR_TEST_H
#define BOND_DISTRIBUTOR_TEST_H

#include <ddMd/communicate/AtomDistributor.h>
#include <ddMd/communicate/BondDistributor.h>
#include <ddMd/communicate/Domain.h>
#include <ddMd/communicate/Buffer.h>
#include <ddMd/storage/AtomStorage.h>
#include <ddMd/storage/AtomIterator.h>
#include <ddMd/chemistry/Atom.h>
#include <ddMd/storage/BondStorage.h>
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

class BondDistributorTest: public ParamFileTest<BondDistributor>
{
public:

   virtual void setUp()
   {}

   void testDistribute()
   {
      printMethod(TEST_FUNC);

      Boundary boundary;
      Domain   domain;
      Buffer   buffer;
      AtomStorage  atomStorage;
      BondStorage  bondStorage;
      AtomDistributor  atomDistributor;
      BondDistributor  bondDistributor;
      std::ifstream configFile;

      // Set connections between atomDistributors
      domain.setBoundary(boundary);
      atomDistributor.associate(domain, boundary, buffer);
      bondDistributor.associate(domain, atomStorage, bondStorage, buffer);

      #ifdef UTIL_MPI
      // Set communicators
      domain.setGridCommunicator(communicator());
      domain.setParamCommunicator(communicator());
      atomStorage.setParamCommunicator(communicator());
      bondStorage.setParamCommunicator(communicator());
      buffer.setParamCommunicator(communicator());
      atomDistributor.setParamCommunicator(communicator());
      bondDistributor.setParamCommunicator(communicator());
      #else
      domain.setRank(0);
      #endif

      // Open parameter file
      #ifdef UTIL_MPI
      openFile("in/BondDistributor.213");
      #else
      openFile("in/BondDistributor.111");
      #endif

      domain.readParam(file());

      // Allocate memory for received atoms on the local processor
      // Maximum total number of atoms received on one processor
      // int atomCapacity = 100; 
      // atomStorage.setParam(atomCapacity, atomCapacity, 200);
      atomStorage.readParam(file());
      bondStorage.readParam(file());

      // Initialize buffer
      //int sendsize = 10;      
      //buffer.allocate(sendsize , sendsize);
      buffer.readParam(file());

      // Initialize AtomDistributor atomDistributor
      // int cacheCapacity = 35; 
      // atomDistributor.setParam(cacheCapacity);
      atomDistributor.readParam(file());
      bondDistributor.readParam(file());

      // Finish reading parameter file
      closeFile();

      Vector boundarylength(6.0, 3.0, 9.0);
      boundary.setLengths(boundarylength);

      int atomCount = 0; // Number of atoms to be distributed by master
      int bondCount = 0; // Number to bonds be distributed by master
      int myRank    = domain.gridRank();

      // If I am the master processor.
      if (myRank == 0) {
         int i;
         Vector  pos;
         Atom*   ptr;

         configFile.open("in/config");
         // Read Max number of atoms to be distributed by the master processor
         configFile >> atomCount;

         std::cout << std::endl;
         std::cout << "Num Atoms to be distributed = " << atomCount << std::endl;

         #if UTIL_MPI
         // Initialize the sendbuffer.
         atomDistributor.initSendBuffer();
         #endif

         // Fill the atom atomDistributors
         for(i = 0; i < atomCount; ++i) {
            ptr = atomDistributor.newAtomPtr();
            ptr->setId(i);
            ptr->setTypeId(0);

            //Assign a random position within the boundary for the atom.
            configFile >> ptr->position();

            //Use position vector for velocity for now
            ptr->velocity() = ptr->position();

            atomDistributor.addAtom(atomStorage);

         }
         file().close();

         // Send any atoms not sent previously.
         atomDistributor.send();

      } else { // If I am not the master processor

         atomDistributor.receive(atomStorage);

      }

      int recvCount = atomStorage.nAtom();
      AtomIterator iter;
      atomStorage.begin(iter);
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

      // Read bonds
      if (myRank == 0) {

         int i;
         Group<2>* ptr;
         //Group<2>  bond;
         //ptr = &bond;

         // Read Max number of atoms to be distributed by the master processor
         configFile >> bondCount;

         std::cout << std::endl;
         std::cout << "Num Bonds to be distributed = " << bondCount << std::endl;

         #if UTIL_MPI
         // Initialize the sendbuffer.
         bondDistributor.initSendBuffer();
         #endif

         // Fill the atom atomDistributors
         for(i = 0; i < bondCount; ++i) {
            ptr = bondDistributor.newPtr();

            //Assign a random position within the boundary for the atom.
            configFile >> *ptr;

            #if 0
            std::cout << "   id "     << ptr->id()
                      << ",  typeId " << ptr->typeId() 
                      << ",  atom0 " << ptr->atomId(0) 
                      << ",  atom1 " << ptr->atomId(1) 
                      << std::endl;
            #endif
       
            bondDistributor.add();

         }
         file().close();

         // Send any bonds not sent previously.
         bondDistributor.send();

      } else { // If I am not the master processor

         bondDistributor.receive();

      }

   }
};

TEST_BEGIN(BondDistributorTest)
TEST_ADD(BondDistributorTest, testDistribute)
TEST_END(BondDistributorTest)

#endif /* BOND_DISTRIBUTOR_TEST_H */
