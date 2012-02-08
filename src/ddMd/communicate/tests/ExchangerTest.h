#ifndef EXCHANGER_TEST_H
#define EXCHANGER_TEST_H

#ifdef UTIL_MPI

#include <ddMd/communicate/Exchanger.h>
#include <ddMd/communicate/AtomDistributor.h>
#include <ddMd/communicate/BondDistributor.h>
#include <ddMd/communicate/Domain.h>
#include <ddMd/communicate/Buffer.h>
#include <ddMd/storage/AtomStorage.h>
#include <ddMd/storage/AtomIterator.h>
#include <ddMd/storage/GhostIterator.h>
#include <ddMd/storage/BondStorage.h>
#include <ddMd/chemistry/Atom.h>
#include <util/space/Grid.h>
#include <util/space/IntVector.h>
#include <util/random/Random.h>
#include <util/mpi/MpiLogger.h>

#ifndef TEST_MPI
#define TEST_MPI
#endif

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>
#include <test/ParamFileTest.h>

using namespace Util;
using namespace DdMd;

class ExchangerTest: public ParamFileTest<Exchanger>
{

      Boundary boundary;
      Domain domain;
      Buffer buffer;
      AtomStorage atomStorage;
      BondStorage bondStorage;
      AtomDistributor atomDistributor;
      BondDistributor bondDistributor;
      Random random;
      int atomCount;

public:

   virtual void setUp()
   {

      // Set connections between objects
      domain.setBoundary(boundary);
      atomDistributor.associate(boundary, domain, buffer);
      bondDistributor.associate(bondStorage, atomStorage, domain, buffer);
      object().associate(boundary, domain, 
                         atomStorage, bondStorage, buffer);

      // Set communicators
      domain.setGridCommunicator(communicator());
      domain.setParamCommunicator(communicator());
      buffer.setParamCommunicator(communicator());
      atomStorage.setParamCommunicator(communicator());
      bondStorage.setParamCommunicator(communicator());
      atomDistributor.setParamCommunicator(communicator());
      bondDistributor.setParamCommunicator(communicator());
      random.setParamCommunicator(communicator());

      // Open parameter file
      openFile("in/Exchanger");

      // Initialize Domain
      domain.readParam(file());

      // Initialize AtomStorage
      atomStorage.readParam(file());
      bondStorage.readParam(file());

      // Initialize Buffer
      buffer.readParam(file());

      // Initialize AtomDistributor 
      atomDistributor.readParam(file());
      bondDistributor.readParam(file());

      // Initialize Random seed
      random.readParam(file());

      object().allocate();

      // Finish reading parameter file
      closeFile();

      std::ifstream configFile;
      int myRank = domain.gridRank();

      if (myRank == 0) {
         configFile.open("in/config");
      }

      // Read and broadcast boundary
      if (myRank == 0) {

         // Read and broadcast system Boundary 
         configFile >> Label("BOUNDARY");
         configFile >> boundary;
      } 
      #if UTIL_MPI
      // Receive broadcast of boundary
      bcast(domain.communicator(), boundary, 0);
      #endif

      // Read and distribute atoms 
      if (myRank == 0) {

         // Read and distribute atoms
         configFile >> Label("ATOMS");

         // Read number of atoms
         int nAtom; 
         configFile >> Label("nAtom") >> nAtom;

         std::cout << std::endl;
         std::cout << "Num Atoms to be distributed = " 
                   << nAtom << std::endl;

         #if UTIL_MPI
         //Initialize the send buffer.
         atomDistributor.initSendBuffer();
         #endif

         // Read atoms
         Atom* atomPtr;
         int id, typeId;
         for(int i = 0; i < nAtom; ++i) {

            // Get pointer to new atom in distributor memory.
            atomPtr = atomDistributor.newAtomPtr();

            configFile >> id >> typeId;
            atomPtr->setId(id);
            atomPtr->setTypeId(typeId);
            configFile >> atomPtr->position();

            // Use position vector for velocity for now
            atomPtr->velocity() = atomPtr->position();

            // Add atom to list for sending.
            atomDistributor.addAtom(atomStorage);

         }

         // Send any atoms not sent previously.
         atomDistributor.send();

      } else { // If I am not the master processor

         #if UTIL_MPI
         // Receive all atoms into AtomStorage
         atomDistributor.receive(atomStorage);
         #endif

      }

      // Check that all atoms are accounted for after distribution.
      {
         int nAtomLocal = atomStorage.nAtom();
         int nAtomAll;
         #ifdef UTIL_MPI
         domain.communicator().Reduce(&nAtomLocal, &nAtomAll, 1, 
                                        MPI::INT, MPI::SUM, 0);
         #else
         nAtomAll = nAtomLocal;
         #endif
         if (myRank == 0) {
            if (nAtomAll != nAtom) {
               UTIL_THROW("nAtomAll != nAtom after distribution");
            }
            std::cout << "nAtomAll after distribution " << nAtomAll 
                      << std::endl;
         }
      }

      // Read and distribute bonds
      if (myRank == 0) {

         // Read and distribute bonds
         configFile >> Label("BONDS");

         // Read number of bonds
         int nBond;  
         configFile >> Label("nBond") >> nBond;

         std::cout << std::endl;
         std::cout << "Num Bonds to be distributed = " 
                   << nBond << std::endl;

         #if UTIL_MPI
         //Initialize the send buffer.
         bondDistributor.initSendBuffer();
         #endif

         // Fill the bond objects
         Bond* bondPtr;
         for (int i = 0; i < nBond; ++i) {

            bondPtr = bondDistributor.newPtr();
            configFile >> *bondPtr;
            bondDistributor.add();

         }

         // Send any bonds not sent previously.
         bondDistributor.send();
      
      } else { // If I am not the master processor

         // Receive all bonds into BondStorage
         bondDistributor.receive();

      }

      if (myRank == 0) {
         configFile.close();
      }

      #if 0
      // Begin distributing atoms 
      Vector boundarylength(6.0, 3.0, 9.0);
      boundary.setLengths(boundarylength);

      atomCount = 0;  // Number to be distributed by master
      int myRank  = domain.gridRank();

      // If I am the master processor.
      if (myRank == 0) {

         std::ifstream atomPosFile;
         Vector pos;
         Atom *ptr;
         int i;

         atomPosFile.open("in/Atompositions");
         // Read Max number of atoms to be distributed by the master processor
         atomPosFile >> atomCount;

         std::cout << std::endl;
         std::cout << "Num Atoms to be distributed = " << atomCount
                   << std::endl;

         // Initialize the send buffer.
         atomDistributor.initSendBuffer();

         // Fill the atom objects
         for(i = 0; i < atomCount; ++i) {

            ptr = atomDistributor.newAtomPtr();

            ptr->setId(i);
            ptr->setTypeId(0);
            atomPosFile >> ptr->position();

            // Use position vector for velocity for now
            ptr->velocity() = ptr->position();

            atomDistributor.addAtom(atomStorage);

         }
         file().close();

         // Send any atoms not sent previously.
         atomDistributor.send();

      } else { // If I am not the master processor

         atomDistributor.receive(atomStorage);

      }
      #endif

   }

   void exchangeAtoms()
   {
      // Range of random increments for the atom positions.
      double range1 = -0.3;
      double range2 = 0.3;

      // Add random increments to atom positions
      AtomIterator  atomIter;
      for(int i = 0; i < 3; i++) {
         atomStorage.begin(atomIter);
         for ( ; !atomIter.atEnd(); ++atomIter) {
            atomIter->position()[i] += random.uniform(range1, range2);
         }
      }

      // Exchange atoms among processors.
      object().exchangeAtoms();

   }

   void testAtomExchange()
   {
      printMethod(TEST_FUNC);

      int  nAtom = 0;     // Number received on this processor.
      int  nAtomAll  = 0; // Number received on all processors.
      int  myRank = domain.gridRank();

      // Check that all atoms are accounted for after distribution.
      nAtom = atomStorage.nAtom();
      communicator().Reduce(&nAtom, &nAtomAll, 1, MPI::INT, MPI::SUM, 0);
      if (myRank == 0) {
         std::cout << std::endl;
         //std::cout << "Total atom count (post-distribute) = " 
         //          << nAtomAll << std::endl;
         TEST_ASSERT(nAtomAll == atomCount);
      }

      // Check that all atoms are within the processor domain.
      AtomIterator  atomIter;
      atomStorage.begin(atomIter);
      for ( ; !atomIter.atEnd(); ++atomIter) {
         TEST_ASSERT(domain.isInDomain(atomIter->position()));
      }

      MpiLogger logger;

      #if 0
      logger.begin();
      std::cout << "Processor: " << myRank << ", Post-distribute nAtom = "
                << atomStorage.nAtom() << std::endl;
      logger.end();
      #endif

      exchangeAtoms();

      // Check that all atoms are accounted for after exchange.
      nAtom = atomStorage.nAtom();
      communicator().Reduce(&nAtom, &nAtomAll, 1, MPI::INT, MPI::SUM, 0);
      if (myRank == 0) {
         std::cout << "Total atom count (post atom exchange) = " << nAtomAll << std::endl;
         TEST_ASSERT(nAtomAll == atomCount);
      }

      // Check that all atoms are within the processor domain.
      atomStorage.begin(atomIter);
      for ( ; !atomIter.atEnd(); ++atomIter) {
         TEST_ASSERT(domain.isInDomain(atomIter->position()));
      }

      TEST_ASSERT(atomStorage.isValid());

      #if 0
      //Print the number of atoms with each processor after the exchange.
      logger.begin();
      std::cout << "Processor " << myRank << " : Post-exchange Atoms count = "
                << atomStorage.nAtom() << std::endl;
      logger.end();
      #endif

   }

   void testGhostExchange()
   {
      printMethod(TEST_FUNC);

      int  nAtom = 0;     // Number received on this processor.
      int  nAtomAll  = 0; // Number received on all processors.
      int  myRank = domain.gridRank();

      AtomIterator  atomIter;
      GhostIterator ghostIter;

      exchangeAtoms();

      // Record number of atoms before exchange
      nAtom = atomStorage.nAtom();

      // Setup ghost exchange
      double pairCutoff = double(0.3);
      object().setPairCutoff(pairCutoff);

      // Exchange ghosts among processors.
      // Vector length = boundary.lengths();
      object().exchangeGhosts();

      // Check that the number of atoms on each processor is unchanged.
      TEST_ASSERT(nAtom == atomStorage.nAtom());

      // Check that all atoms are accounted for after atom and ghost exchanges.
      nAtom = atomStorage.nAtom();
      communicator().Reduce(&nAtom, &nAtomAll, 1, MPI::INT, MPI::SUM, 0);
      if (myRank == 0) {
         // std::cout << "Total atom count (post ghost exchange) = " 
         //           << nAtomAll << std::endl;
         TEST_ASSERT(nAtomAll == atomCount);
      }

      // Check that all local atoms are within the processor domain.
      atomStorage.begin(atomIter);
      for ( ; !atomIter.atEnd(); ++atomIter) {
         TEST_ASSERT(domain.isInDomain(atomIter->position()));
      }

      // Check that all ghosts are outside the processor domain.
      atomStorage.begin(ghostIter);
      for ( ; !ghostIter.atEnd(); ++ghostIter) {
         TEST_ASSERT(!domain.isInDomain(ghostIter->position()));
      }

      TEST_ASSERT(atomStorage.isValid());

      #if 0
      MpiLogger logger;

      //Print number of atoms on each processor after the ghost exchange.
      logger.begin();
      std::cout << "Processor " << myRank 
                << " : Post-ghost exchange Atom  count = "
                << atomStorage.nAtom() << std::endl;
      logger.end();

      // Print number of ghosts on each processor after the exchange.
      logger.begin();
      std::cout << "Processor " << myRank 
                << " : Post-ghost exchange Ghost count = "
                << atomStorage.nGhost() << std::endl;
      logger.end();
      #endif

   }

   void testGhostUpdate()
   {
      printMethod(TEST_FUNC);

      int  nAtom  = 0;    // Number of atoms on this processor.
      int  nGhost = 0;    // Number of ghosts on this processor.
      int  nAtomAll  = 0; // Number received on all processors.
      int  myRank = domain.gridRank();

      AtomIterator   atomIter;
      GhostIterator  ghostIter;
      DArray<Vector> ghostPositions;

      exchangeAtoms();

      // Record number of atoms before ghost exchange
      nAtom = atomStorage.nAtom();

      // Set slab width used for ghost exchange.
      double pairCutoff = double(0.3);
      object().setPairCutoff(pairCutoff);

      // Exchange ghosts among processors.
      // Vector length = boundary.lengths();
      object().exchangeGhosts();

      // Record number of ghosts after ghost exchange, before update.
      nGhost = atomStorage.nGhost();

      // Update ghost positions
      object().updateGhosts();

      // Check that the number of atoms on each processor is unchanged.
      TEST_ASSERT(nAtom == atomStorage.nAtom());

      // Check that the number of ghosts on each processor is unchanged.
      TEST_ASSERT(nGhost == atomStorage.nGhost());

      // Check that all atoms are accounted for after atom and ghost exchanges.
      nAtom = atomStorage.nAtom();
      communicator().Reduce(&nAtom, &nAtomAll, 1, MPI::INT, MPI::SUM, 0);
      if (myRank == 0) {
         // std::cout << "Total atom count (post ghost exchange) = " 
         //           << nAtomAll << std::endl;
         TEST_ASSERT(nAtomAll == atomCount);
      }

      // Check that all atoms are within the processor domain.
      atomStorage.begin(atomIter);
      for ( ; !atomIter.atEnd(); ++atomIter) {
         TEST_ASSERT(domain.isInDomain(atomIter->position()));
      }

      // Check that all ghosts are outside the processor domain.
      atomStorage.begin(ghostIter);
      for ( ; !ghostIter.atEnd(); ++ghostIter) {
         TEST_ASSERT(!domain.isInDomain(ghostIter->position()));
      }

      TEST_ASSERT(atomStorage.isValid());

      #if 0
      MpiLogger logger;

      //Print number of atoms on each processor after the ghost exchange.
      logger.begin();
      std::cout << "Processor " << myRank << " : Post-ghost exchange Atom  count = "
                << atomStorage.nAtom() << std::endl;
      logger.end();

      // Print number of ghosts on each processor after the exchange.
      logger.begin();
      std::cout << "Processor " << myRank << " : Post-ghost exchange Ghost count = "
                << atomStorage.nGhost() << std::endl;
      logger.end();
      #endif

   }

   void testGhostUpdateCycle()
   {
      printMethod(TEST_FUNC);

      int  nAtom  = 0;    // Number of atoms on this processor.
      int  nGhost = 0;    // Number of ghosts on this processor.
      int  nAtomAll  = 0; // Number received on all processors.
      int  myRank = domain.gridRank();

      AtomIterator   atomIter;
      GhostIterator  ghostIter;
      DArray<Vector> ghostPositions;

      // Set slab width used for ghost exchange.
      double pairCutoff = double(0.5);
      object().setPairCutoff(pairCutoff);

      for (int i=0; i < 3; ++i) {

         // Clear ghosts prior to exchanging atoms.
         // atomStorage.clearGhosts();

         // Move positions and exchange ownership
         exchangeAtoms();

         TEST_ASSERT(atomStorage.isValid());

         // Record number of atoms before ghost exchange
         nAtom = atomStorage.nAtom();
   
         // Exchange ghosts among processors.
         // Vector length = boundary.lengths();
         object().exchangeGhosts();
   
         TEST_ASSERT(atomStorage.isValid());

         // Record number of ghosts after ghost exchange, before update.
         nGhost = atomStorage.nGhost();
   
         for (int j=0; j < 3; ++j) {

            // Update ghost positions
            object().updateGhosts();

            TEST_ASSERT(atomStorage.isValid());
      
            // Assert number of atoms on each processor is unchanged.
            TEST_ASSERT(nAtom == atomStorage.nAtom());
      
            // Assert number of ghosts on each processor is unchanged.
            TEST_ASSERT(nGhost == atomStorage.nGhost());
      
            // Assert atoms accounted for after atom & ghost exchanges.
            nAtom = atomStorage.nAtom();
            communicator().Reduce(&nAtom, &nAtomAll, 1, MPI::INT, MPI::SUM, 0);
            if (myRank == 0) {
               // std::cout << "Total atom count (post ghost exchange) = " 
               //           << nAtomAll << std::endl;
               TEST_ASSERT(nAtomAll == atomCount);
            }
      
            // Assert all atoms are within the processor domain.
            atomStorage.begin(atomIter);
            for ( ; !atomIter.atEnd(); ++atomIter) {
               TEST_ASSERT(domain.isInDomain(atomIter->position()));
            }
      
            // Assert all ghosts are outside the processor domain.
            atomStorage.begin(ghostIter);
            for ( ; !ghostIter.atEnd(); ++ghostIter) {
               TEST_ASSERT(!domain.isInDomain(ghostIter->position()));
            }
   
            #if 0
            MpiLogger logger;
      
            //Print number atoms on each node after ghost exchange.
            logger.begin();
            std::cout << "Processor " << myRank 
                      << " : Post-ghost exchange Atom  count = "
                      << atomStorage.nAtom() << std::endl;
            logger.end();
      
            // Print number ghosts on each node after the exchange.
            logger.begin();
            std::cout << "Processor " << myRank 
                      << " : Post-ghost exchange Ghost count = "
                      << atomStorage.nGhost() << std::endl;
            logger.end();
            #endif

         } // end ghost cycle update

      } // end exchange update

   }

};

TEST_BEGIN(ExchangerTest)
TEST_ADD(ExchangerTest, testAtomExchange)
TEST_ADD(ExchangerTest, testGhostExchange)
TEST_ADD(ExchangerTest, testGhostUpdate)
TEST_ADD(ExchangerTest, testGhostUpdateCycle)
TEST_END(ExchangerTest)

#endif /* UTIL_MPI */
#endif /* EXCHANGER_TEST_H */
