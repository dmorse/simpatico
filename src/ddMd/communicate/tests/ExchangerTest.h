#ifndef EXCHANGER_TEST_H
#define EXCHANGER_TEST_H

#ifdef UTIL_MPI

#include <ddMd/communicate/Exchanger.h>
#include <ddMd/communicate/AtomDistributor.h>
#include <ddMd/communicate/Domain.h>
#include <ddMd/communicate/Buffer.h>
#include <ddMd/storage/AtomStorage.h>
#include <ddMd/storage/AtomIterator.h>
#include <ddMd/storage/GhostIterator.h>
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

      Boundary    boundary;
      Domain      domain;
      AtomStorage     storage;
      Buffer      buffer;
      AtomDistributor distributor;
      Random      random;
      int         atomCount;

public:

   virtual void setUp()
   {

      // Set connections between objects
      domain.setBoundary(boundary);
      distributor.associate(boundary, domain, buffer);
      object().associate(boundary, domain, storage, buffer);

      // Set communicators
      domain.setGridCommunicator(communicator());
      domain.setParamCommunicator(communicator());
      buffer.setParamCommunicator(communicator());
      storage.setParamCommunicator(communicator());
      distributor.setParamCommunicator(communicator());
      random.setParamCommunicator(communicator());

      // Open parameter file
      openFile("in/Exchanger");

      // Initialize Domain
      domain.readParam(file());

      // Initialize AtomStorage
      storage.readParam(file());

      // Initialize Buffer
      buffer.readParam(file());

      // Initialize AtomDistributor 
      distributor.readParam(file());

      // Initialize Random seed
      random.readParam(file());

      object().allocate();

      // Finish reading parameter file
      closeFile();

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
         distributor.initSendBuffer();

         // Fill the atom objects
         for(i = 0; i < atomCount; ++i) {

            ptr = distributor.newAtomPtr();

            ptr->setId(i);
            ptr->setTypeId(0);
            atomPosFile >> ptr->position();

            // Use position vector for velocity for now
            ptr->velocity() = ptr->position();

            distributor.addAtom(storage);

         }
         file().close();

         // Send any atoms not sent previously.
         distributor.send();

      } else { // If I am not the master processor

         distributor.receive(storage);

      }

   }

   void exchangeAtoms()
   {
      // Range of random increments for the atom positions.
      double range1 = double(-0.3);
      double range2 = double(0.3);

      // Add random increments to atom positions
      AtomIterator  atomIter;
      for(int i = 0; i < 3; i++) {
        storage.begin(atomIter);
    	 for ( ; !atomIter.atEnd(); ++atomIter) {
    	    atomIter->position()[i] += random.uniform(range1, range2);
     	 }
      }

      // Exchange atoms among processors.
      // Vector length = boundary.lengths();
      object().exchangeAtoms();

   }

   void testAtomExchange()
   {
      printMethod(TEST_FUNC);

      int  nAtom = 0;     // Number received on this processor.
      int  nAtomAll  = 0; // Number received on all processors.
      int  myRank = domain.gridRank();


      // Check that all atoms are accounted for after distribution.
      nAtom = storage.nAtom();
      communicator().Reduce(&nAtom, &nAtomAll, 1, MPI::INT, MPI::SUM, 0);
      if (myRank == 0) {
         std::cout << std::endl;
         //std::cout << "Total atom count (post-distribute) = " 
         //          << nAtomAll << std::endl;
         TEST_ASSERT(nAtomAll == atomCount);
      }

      // Check that all atoms are within the processor domain.
      AtomIterator  atomIter;
      storage.begin(atomIter);
      for ( ; !atomIter.atEnd(); ++atomIter) {
         TEST_ASSERT(domain.isInDomain(atomIter->position()));
      }

      MpiLogger logger;

      #if 0
      logger.begin();
      std::cout << "Processor: " << myRank << ", Post-distribute nAtom = "
                << storage.nAtom() << std::endl;
      logger.end();
      #endif

      exchangeAtoms();

      // Check that all atoms are accounted for after exchange.
      nAtom = storage.nAtom();
      communicator().Reduce(&nAtom, &nAtomAll, 1, MPI::INT, MPI::SUM, 0);
      if (myRank == 0) {
         std::cout << "Total atom count (post atom exchange) = " << nAtomAll << std::endl;
         TEST_ASSERT(nAtomAll == atomCount);
      }

      // Check that all atoms are within the processor domain.
      storage.begin(atomIter);
      for ( ; !atomIter.atEnd(); ++atomIter) {
         TEST_ASSERT(domain.isInDomain(atomIter->position()));
      }

      TEST_ASSERT(storage.isValid());

      #if 0
      //Print the number of atoms with each processor after the exchange.
      logger.begin();
      std::cout << "Processor " << myRank << " : Post-exchange Atoms count = "
                << storage.nAtom() << std::endl;
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
      nAtom = storage.nAtom();

      // Setup ghost exchange
      double pairCutoff = double(0.3);
      object().setPairCutoff(pairCutoff);

      // Exchange ghosts among processors.
      // Vector length = boundary.lengths();
      object().exchangeGhosts();

      // Check that the number of atoms on each processor is unchanged.
      TEST_ASSERT(nAtom == storage.nAtom());

      // Check that all atoms are accounted for after atom and ghost exchanges.
      nAtom = storage.nAtom();
      communicator().Reduce(&nAtom, &nAtomAll, 1, MPI::INT, MPI::SUM, 0);
      if (myRank == 0) {
         // std::cout << "Total atom count (post ghost exchange) = " 
         //           << nAtomAll << std::endl;
         TEST_ASSERT(nAtomAll == atomCount);
      }

      // Check that all local atoms are within the processor domain.
      storage.begin(atomIter);
      for ( ; !atomIter.atEnd(); ++atomIter) {
         TEST_ASSERT(domain.isInDomain(atomIter->position()));
      }

      // Check that all ghosts are outside the processor domain.
      storage.begin(ghostIter);
      for ( ; !ghostIter.atEnd(); ++ghostIter) {
         TEST_ASSERT(!domain.isInDomain(ghostIter->position()));
      }

      TEST_ASSERT(storage.isValid());

      #if 0
      MpiLogger logger;

      //Print number of atoms on each processor after the ghost exchange.
      logger.begin();
      std::cout << "Processor " << myRank 
                << " : Post-ghost exchange Atom  count = "
                << storage.nAtom() << std::endl;
      logger.end();

      // Print number of ghosts on each processor after the exchange.
      logger.begin();
      std::cout << "Processor " << myRank 
                << " : Post-ghost exchange Ghost count = "
                << storage.nGhost() << std::endl;
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
      nAtom = storage.nAtom();

      // Set slab width used for ghost exchange.
      double pairCutoff = double(0.3);
      object().setPairCutoff(pairCutoff);

      // Exchange ghosts among processors.
      // Vector length = boundary.lengths();
      object().exchangeGhosts();

      // Record number of ghosts after ghost exchange, before update.
      nGhost = storage.nGhost();

      // Update ghost positions
      object().updateGhosts();

      // Check that the number of atoms on each processor is unchanged.
      TEST_ASSERT(nAtom == storage.nAtom());

      // Check that the number of ghosts on each processor is unchanged.
      TEST_ASSERT(nGhost == storage.nGhost());

      // Check that all atoms are accounted for after atom and ghost exchanges.
      nAtom = storage.nAtom();
      communicator().Reduce(&nAtom, &nAtomAll, 1, MPI::INT, MPI::SUM, 0);
      if (myRank == 0) {
         // std::cout << "Total atom count (post ghost exchange) = " 
         //           << nAtomAll << std::endl;
         TEST_ASSERT(nAtomAll == atomCount);
      }

      // Check that all atoms are within the processor domain.
      storage.begin(atomIter);
      for ( ; !atomIter.atEnd(); ++atomIter) {
         TEST_ASSERT(domain.isInDomain(atomIter->position()));
      }

      // Check that all ghosts are outside the processor domain.
      storage.begin(ghostIter);
      for ( ; !ghostIter.atEnd(); ++ghostIter) {
         TEST_ASSERT(!domain.isInDomain(ghostIter->position()));
      }

      TEST_ASSERT(storage.isValid());

      #if 0
      MpiLogger logger;

      //Print number of atoms on each processor after the ghost exchange.
      logger.begin();
      std::cout << "Processor " << myRank << " : Post-ghost exchange Atom  count = "
                << storage.nAtom() << std::endl;
      logger.end();

      // Print number of ghosts on each processor after the exchange.
      logger.begin();
      std::cout << "Processor " << myRank << " : Post-ghost exchange Ghost count = "
                << storage.nGhost() << std::endl;
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
         storage.clearGhosts();

         // Move positions and exchange ownership
         exchangeAtoms();

         TEST_ASSERT(storage.isValid());

         // Record number of atoms before ghost exchange
         nAtom = storage.nAtom();
   
         // Exchange ghosts among processors.
         // Vector length = boundary.lengths();
         object().exchangeGhosts();
   
         TEST_ASSERT(storage.isValid());

         // Record number of ghosts after ghost exchange, before update.
         nGhost = storage.nGhost();
   
         for (int j=0; j < 3; ++j) {

            // Update ghost positions
            object().updateGhosts();

            TEST_ASSERT(storage.isValid());
      
            // Assert number of atoms on each processor is unchanged.
            TEST_ASSERT(nAtom == storage.nAtom());
      
            // Assert number of ghosts on each processor is unchanged.
            TEST_ASSERT(nGhost == storage.nGhost());
      
            // Assert atoms accounted for after atom & ghost exchanges.
            nAtom = storage.nAtom();
            communicator().Reduce(&nAtom, &nAtomAll, 1, MPI::INT, MPI::SUM, 0);
            if (myRank == 0) {
               // std::cout << "Total atom count (post ghost exchange) = " 
               //           << nAtomAll << std::endl;
               TEST_ASSERT(nAtomAll == atomCount);
            }
      
            // Assert all atoms are within the processor domain.
            storage.begin(atomIter);
            for ( ; !atomIter.atEnd(); ++atomIter) {
               TEST_ASSERT(domain.isInDomain(atomIter->position()));
            }
      
            // Assert all ghosts are outside the processor domain.
            storage.begin(ghostIter);
            for ( ; !ghostIter.atEnd(); ++ghostIter) {
               TEST_ASSERT(!domain.isInDomain(ghostIter->position()));
            }
   
            #if 0
            MpiLogger logger;
      
            //Print number atoms on each node after ghost exchange.
            logger.begin();
            std::cout << "Processor " << myRank 
                      << " : Post-ghost exchange Atom  count = "
                      << storage.nAtom() << std::endl;
            logger.end();
      
            // Print number ghosts on each node after the exchange.
            logger.begin();
            std::cout << "Processor " << myRank 
                      << " : Post-ghost exchange Ghost count = "
                      << storage.nGhost() << std::endl;
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
