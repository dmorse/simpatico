#ifndef DDMD_DISTRIBUTOR_TEST_H
#define DDMD_DISTRIBUTOR_TEST_H

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

class AtomDistributorTest: public ParamFileTest
{
private:
   
   AtomDistributor distributor_;
   Boundary boundary_;
   Domain domain_;
   Buffer buffer_;
   AtomStorage storage_;

public:

   virtual void setUp()
   {

      // Set connections between objects
      domain_.setBoundary(boundary_);
      distributor_.associate(domain_, boundary_, storage_, buffer_);

      #ifdef UTIL_MPI
      // Set communicators
      domain_.setGridCommunicator(communicator());
      domain_.setParamCommunicator(communicator());
      storage_.setParamCommunicator(communicator());
      buffer_.setParamCommunicator(communicator());
      distributor_.setParamCommunicator(communicator());
      #else
      domain_.setRank(0);
      #endif

   }

   void testDistribute()
   {
      printMethod(TEST_FUNC);
      std::ifstream atomposfile;

      // Open parameter file
      #ifdef UTIL_MPI
      openFile("in/AtomDistributor.213");
      #else
      openFile("in/AtomDistributor.111");
      #endif

      domain_.readParam(file());

      // Allocate memory for received atoms on the local processor
      // Maximum total number of atoms received on one processor
      // int atomCapacity = 100; 
      // storage_.initialize(atomCapacity, atomCapacity, 200);
      storage_.readParam(file());

      // Initialize buffer
      //int sendsize = 10;      
      //buffer_.allocate(sendsize , sendsize);
      buffer_.readParam(file());

      // Initialize AtomDistributor object
      // int cacheCapacity = 35; 
      // distributor_.initialize(cacheCapacity);
      distributor_.readParam(file());

      // Finish reading parameter file
      closeFile();

      Vector boundarylength(6.0, 3.0, 9.0);
      boundary_.setOrthorhombic(boundarylength);

      int atomCount = 0; // Number to be distributed by master
      int myRank    = domain_.gridRank();

      // If I am the master processor.
      if (myRank == 0) {
         int i;
         Vector  pos;
         Atom*   ptr;

         openInputFile("in/Atompositions", atomposfile);
         //atomposfile.open("in/Atompositions");
         // Read Max number of atoms to be distributed by the master processor
         atomposfile >> atomCount;

         //std::cout << std::endl;
         //std::cout << "Num Atoms to be distributed = " 
         //          << atomCount << std::endl;

         #if UTIL_MPI
         // Initialize the sendbuffer.
         distributor_.setup();
         #endif

         // Fill the atom objects
         for(i = 0; i < atomCount; ++i) {
            ptr = distributor_.newAtomPtr();
            ptr->setId(i);
            ptr->setTypeId(0);

            // Read a position from file.
            if (UTIL_ORTHOGONAL) {
               atomposfile >> ptr->position();
            } else {
               Vector r;
               atomposfile >> r;
               boundary_.transformCartToGen(r, ptr->position());
            }

            //Use position vector for velocity for now
            ptr->velocity() = ptr->position();

            distributor_.addAtom();

         }
         file().close();

         // Send any atoms not sent previously.
         distributor_.send();

      } else { // If I am not the master processor

         distributor_.receive();

      }

      int recvCount = storage_.nAtom();
      AtomIterator iter;
      storage_.begin(iter);
      for ( ; iter.notEnd(); ++iter) {
         TEST_ASSERT(domain_.isInDomain(iter->position()));
      }

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
};

TEST_BEGIN(AtomDistributorTest)
TEST_ADD(AtomDistributorTest, testDistribute)
TEST_END(AtomDistributorTest)

#endif /* DISTRIBUTOR_TEST_H */
