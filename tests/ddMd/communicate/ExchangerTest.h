#ifndef DDMD_EXCHANGER_TEST_H
#define DDMD_EXCHANGER_TEST_H

#include <ddMd/configIos/DdMdConfigIo.h>
#include <ddMd/communicate/Domain.h>
#include <ddMd/communicate/Buffer.h>
#include <ddMd/communicate/Exchanger.h>
#include <ddMd/storage/AtomStorage.h>
#include <ddMd/storage/AtomIterator.h>
#include <ddMd/storage/GhostIterator.h>
#include <ddMd/storage/BondStorage.h>
#ifdef INTER_ANGLE
#include <ddMd/storage/AngleStorage.h>
#endif
#ifdef INTER_DIHEDRAL
#include <ddMd/storage/DihedralStorage.h>
#endif
#include <util/boundary/Boundary.h>
#include <ddMd/chemistry/MaskPolicy.h>
#include <util/random/Random.h>
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

class ExchangerTest: public ParamFileTest<Exchanger>
{
private:

   Boundary boundary;
   Domain domain;
   Buffer buffer;
   AtomStorage atomStorage;
   BondStorage bondStorage;
   #ifdef INTER_ANGLE
   AngleStorage angleStorage;
   #endif
   #ifdef INTER_DIHEDRAL
   DihedralStorage dihedralStorage;
   #endif
   DdMdConfigIo configIo;
   Random random;
   int atomCount;

public:

   void setUp()
   {

      // Set connections between atomDistributors
      domain.setBoundary(boundary);
      configIo.associate(domain, boundary, atomStorage, bondStorage, 
                         #ifdef INTER_ANGLE
                         angleStorage,
                         #endif
                         #ifdef INTER_DIHEDRAL
                         dihedralStorage,
                         #endif
                         buffer);
      object().associate(domain, boundary, atomStorage, bondStorage, 
                         #ifdef INTER_ANGLE
                         angleStorage,
                         #endif
                         #ifdef INTER_DIHEDRAL
                         dihedralStorage,
                         #endif
                         buffer);

      #ifdef UTIL_MPI
      // Set communicators
      domain.setGridCommunicator(communicator());
      domain.setParamCommunicator(communicator());
      atomStorage.setParamCommunicator(communicator());
      bondStorage.setParamCommunicator(communicator());
      #ifdef INTER_ANGLE
      angleStorage.setParamCommunicator(communicator());
      #endif
      #ifdef INTER_DIHEDRAL
      dihedralStorage.setParamCommunicator(communicator());
      #endif
      buffer.setParamCommunicator(communicator());
      configIo.setParamCommunicator(communicator());
      random.setParamCommunicator(communicator());
      #else
      domain.setRank(0);
      #endif

      // Open parameter file
      openFile("in/Exchanger");

      domain.readParam(file());
      buffer.readParam(file());
      configIo.readParam(file());
      random.readParam(file());
      atomStorage.readParam(file());
      bondStorage.readParam(file());
      #ifdef INTER_ANGLE
      angleStorage.readParam(file());
      #endif
      #ifdef INTER_DIHEDRAL
      dihedralStorage.readParam(file());
      #endif

      // Finish reading parameter file
      closeFile();

      object().setPairCutoff(0.5);
      object().allocate();

      MaskPolicy policy = MaskBonded;
      std::ifstream configFile("in/config");
      configIo.readConfig(configFile, policy);

      int  nAtom = 0;     // Number received on this processor.
      int  nAtomAll  = 0; // Number received on all processors.

      // Check that all atoms are accounted for after distribution.
      nAtom = atomStorage.nAtom();
      communicator().Reduce(&nAtom, &nAtomAll, 1, MPI::INT, MPI::SUM, 0);
      if (domain.gridRank() == 0) {
         atomCount = nAtomAll;
      }

   }

   void displaceAtoms(double range)
   {
      Vector ranges;

      // Set ranges for random diplacements in each direction.
      if (UTIL_ORTHOGONAL || atomStorage.isCartesian()) {
         for (int i = 0; i < Dimension; ++i) {
            ranges[i] = range;
         }
      } else {
         for (int i = 0; i < Dimension; ++i) {
            ranges[i] = range/boundary.length(i);
         }
      }

      // Iterate over atoms, adding random displacements.
      double min, max;
      AtomIterator atomIter;
      for (atomStorage.begin(atomIter); atomIter.notEnd(); ++atomIter) {
         for (int i = 0; i < Dimension; ++i) {
            max = ranges[i];
            min = -max;
            atomIter->position()[i] += random.uniform(min, max);
         }
      }

   }

   virtual void testDistribute()
   { 
      printMethod(TEST_FUNC); 
   }

   void testAtomExchange()
   {
      printMethod(TEST_FUNC);

      int  nAtom = 0;     // Number received on this processor.
      int  nAtomAll  = 0; // Number received on all processors.
      int  myRank = domain.gridRank();

      // Check that all atoms are within the processor domain.
      AtomIterator  atomIter;
      atomStorage.begin(atomIter);
      for ( ; atomIter.notEnd(); ++atomIter) {
         TEST_ASSERT(domain.isInDomain(atomIter->position()));
      }

      TEST_ASSERT(atomStorage.isValid());
      TEST_ASSERT(bondStorage.isValid(atomStorage, domain.communicator(), 
                  false));

      TEST_ASSERT(!atomStorage.isCartesian());
      double range = 0.4;
      displaceAtoms(range);

      object().exchangeAtoms();

      // Check that all atoms are accounted for after exchange.
      nAtom = atomStorage.nAtom();
      communicator().Reduce(&nAtom, &nAtomAll, 1, MPI::INT, MPI::SUM, 0);
      if (myRank == 0) {
         TEST_ASSERT(nAtomAll == atomCount);
      }

      // Check that all atoms are within the processor domain.
      atomStorage.begin(atomIter);
      for ( ; atomIter.notEnd(); ++atomIter) {
         TEST_ASSERT(domain.isInDomain(atomIter->position()));
      }

      TEST_ASSERT(atomStorage.isValid());
      TEST_ASSERT(bondStorage.isValid(atomStorage, domain.communicator(), 
                  false));

   }

   void testExchange()
   {
      printMethod(TEST_FUNC);

      int  nAtom = 0;     // Number received on this processor.
      int  nAtomAll  = 0; // Number received on all processors.
      int  myRank = domain.gridRank();

      AtomIterator  atomIter;
      GhostIterator ghostIter;

      double range = 0.4;
      displaceAtoms(range);
      
      object().exchange();

      // Check that all atoms are accounted for after ghost exchange.
      nAtom = atomStorage.nAtom();
      communicator().Reduce(&nAtom, &nAtomAll, 1, MPI::INT, MPI::SUM, 0);
      if (myRank == 0) {
         // std::cout << "Total atom count (post ghost exchange) = " 
         //           << nAtomAll << std::endl;
         TEST_ASSERT(nAtomAll == atomCount);
      }

      // Check that all local atoms are within the processor domain.
      atomStorage.begin(atomIter);
      for ( ; atomIter.notEnd(); ++atomIter) {
         TEST_ASSERT(domain.isInDomain(atomIter->position()));
      }

      // Check that all ghosts are outside the processor domain.
      atomStorage.begin(ghostIter);
      for ( ; ghostIter.notEnd(); ++ghostIter) {
         TEST_ASSERT(!domain.isInDomain(ghostIter->position()));
      }

      // Call isVlaid() methods of all storage containers.
      TEST_ASSERT(atomStorage.isValid());
      TEST_ASSERT(bondStorage.isValid(atomStorage, domain.communicator(), 
                  true));
      #ifdef INTER_ANGLE
      TEST_ASSERT(angleStorage.isValid(atomStorage, 
                  domain.communicator(), true));
      #endif
      #ifdef INTER_DIHEDRAL
      TEST_ASSERT(dihedralStorage.isValid(atomStorage, 
                  domain.communicator(), true));
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

      double range = 0.1;
      displaceAtoms(range);

      atomStorage.clearSnapshot();
      object().exchange();

      // Record number of atoms and ghosts after exchange
      nAtom = atomStorage.nAtom();
      nGhost = atomStorage.nGhost();

      // Transform to Cartesian coordinates
      if (!UTIL_ORTHOGONAL) {
         atomStorage.transformGenToCart(boundary);
      }
      atomStorage.makeSnapshot();

      //displaceAtoms(range);

      // Update ghost positions
      object().update();

      // Check number of atoms and ghosts on each processor is unchanged.
      TEST_ASSERT(nAtom == atomStorage.nAtom());
      TEST_ASSERT(nGhost == atomStorage.nGhost());

      // Check that all atoms are accounted for after atom and ghost exchanges.
      communicator().Reduce(&nAtom, &nAtomAll, 1, MPI::INT, MPI::SUM, 0);
      if (myRank == 0) {
         // std::cout << "Total atom count (post ghost exchange) = " 
         //           << nAtomAll << std::endl;
         TEST_ASSERT(nAtomAll == atomCount);
      }

      // Transform back to generalized coordinates
      if (!UTIL_ORTHOGONAL) {
         atomStorage.transformCartToGen(boundary);
      }

      // Check that all atoms are within the processor domain.
      atomStorage.begin(atomIter);
      for ( ; atomIter.notEnd(); ++atomIter) {
         TEST_ASSERT(domain.isInDomain(atomIter->position()));
      }

      // Check that all ghosts are outside the processor domain.
      atomStorage.begin(ghostIter);
      for ( ; ghostIter.notEnd(); ++ghostIter) {
         TEST_ASSERT(!domain.isInDomain(ghostIter->position()));
      }

      TEST_ASSERT(atomStorage.isValid());
      TEST_ASSERT(bondStorage.isValid(atomStorage, domain.communicator(), 
                  true));
      #ifdef INTER_ANGLE
      TEST_ASSERT(angleStorage.isValid(atomStorage, 
                  domain.communicator(), true));
      #endif
      #ifdef INTER_DIHEDRAL
      TEST_ASSERT(dihedralStorage.isValid(atomStorage, 
                  domain.communicator(), true));
      #endif


   }

   void testGhostUpdateCycle()
   {
      printMethod(TEST_FUNC);

      int  nAtom  = 0;    // Number of atoms on this processor.
      int  nGhost = 0;    // Number of ghosts on this processor.

      AtomIterator   atomIter;
      GhostIterator  ghostIter;

      double range = 0.4;
      displaceAtoms(range);

      object().exchange();
      nAtom = atomStorage.nAtom();
      nGhost = atomStorage.nGhost();

      // Assert that all atoms are within the processor domain.
      atomStorage.begin(atomIter);
      for ( ; atomIter.notEnd(); ++atomIter) {
         TEST_ASSERT(domain.isInDomain(atomIter->position()));
      }

      // Assert that all ghosts are outside the processor domain.
      atomStorage.begin(ghostIter);
      for ( ; ghostIter.notEnd(); ++ghostIter) {
         TEST_ASSERT(!domain.isInDomain(ghostIter->position()));
      }

      TEST_ASSERT(atomStorage.isValid());
      TEST_ASSERT(bondStorage.isValid(atomStorage, domain.communicator(), 
                  true));

      // Transform to Cartesian coordinates
      if (!UTIL_ORTHOGONAL) {
         atomStorage.transformGenToCart(boundary);
      }

      range = 0.1;
      for (int i=0; i < 3; ++i) {

         displaceAtoms(range);

         for (int j=0; j < 3; ++j) {
            object().update();
            TEST_ASSERT(nGhost == atomStorage.nGhost());
            TEST_ASSERT(nAtom == atomStorage.nAtom());
            displaceAtoms(range);
         }

         // Transform to Cartesian coordinates
         if (!UTIL_ORTHOGONAL) {
            atomStorage.transformCartToGen(boundary);
         }

         object().exchange();
         nAtom  = atomStorage.nAtom();
         nGhost = atomStorage.nGhost();

         // Check that all atoms are within the processor domain.
         atomStorage.begin(atomIter);
         for ( ; atomIter.notEnd(); ++atomIter) {
            TEST_ASSERT(domain.isInDomain(atomIter->position()));
         }

         // Check that all ghosts are outside the processor domain.
         atomStorage.begin(ghostIter);
         for ( ; ghostIter.notEnd(); ++ghostIter) {
            TEST_ASSERT(!domain.isInDomain(ghostIter->position()));
         }

         TEST_ASSERT(atomStorage.isValid());
         TEST_ASSERT(bondStorage.isValid(atomStorage, domain.communicator(),
                                         true)); 
         #ifdef INTER_ANGLE
         TEST_ASSERT(angleStorage.isValid(atomStorage, 
                     domain.communicator(), true));
         #endif
         #ifdef INTER_DIHEDRAL
         TEST_ASSERT(dihedralStorage.isValid(atomStorage, 
                     domain.communicator(), true));
         #endif

         // Transform to Cartesian coordinates
         if (!UTIL_ORTHOGONAL) {
            atomStorage.transformGenToCart(boundary);
         }

      }

   }

   void testExchangeUpdateCycle()
   {
      printMethod(TEST_FUNC);

      int  nAtom  = 0;    // Number of atoms on this processor.
      int  nGhost = 0;    // Number of ghosts on this processor.

      AtomIterator   atomIter;
      GhostIterator  ghostIter;

      double range = 0.2;

      #if 0
      atomStorage.makeSnapshot();
      displaceAtoms(range);
      #endif

      atomStorage.clearSnapshot();
      object().exchange();

      nAtom = atomStorage.nAtom();
      nGhost = atomStorage.nGhost();

      // Assert that all atoms are within the processor domain.
      atomStorage.begin(atomIter);
      for ( ; atomIter.notEnd(); ++atomIter) {
         TEST_ASSERT(domain.isInDomain(atomIter->position()));
      }

      // Assert that all ghosts are outside the processor domain.
      atomStorage.begin(ghostIter);
      for ( ; ghostIter.notEnd(); ++ghostIter) {
         TEST_ASSERT(!domain.isInDomain(ghostIter->position()));
      }

      // Check Atom and Group Storage containers
      TEST_ASSERT(atomStorage.isValid());
      TEST_ASSERT(bondStorage.isValid(atomStorage, domain.communicator(), 
                  true));
      #ifdef INTER_ANGLE
      TEST_ASSERT(angleStorage.isValid(atomStorage, 
                  domain.communicator(), true));
      #endif
      #ifdef INTER_DIHEDRAL
      TEST_ASSERT(dihedralStorage.isValid(atomStorage, 
                  domain.communicator(), true));
      #endif

      if (!UTIL_ORTHOGONAL) {
         TEST_ASSERT(!atomStorage.isCartesian());
         atomStorage.transformGenToCart(boundary);
      }
      TEST_ASSERT(atomStorage.isCartesian());
      atomStorage.makeSnapshot();

      #if 0
      if (domain.gridRank() == 0) {
         std::cout << std::endl;
      }
      #endif

      range = 0.02;
      double skin = 0.10;
      int  nExchange = 0;
      int  nUpdate = 0;
      int  i, j;
      bool  needExchange;

      j = 0;
      for (i=0; i < 200; ++i) {

         TEST_ASSERT(atomStorage.isCartesian());
         displaceAtoms(range);
         ++j;

         needExchange = atomStorage.needExchange(domain.communicator(), skin);
         if (needExchange || j > 10) {

            #if 0
            if (domain.gridRank() == 0) {
               std::cout << "step i = " << i << "  E" << std::endl;
            }
            #endif

            atomStorage.clearSnapshot();
            if (!UTIL_ORTHOGONAL) {
               atomStorage.transformCartToGen(boundary);
            }
            object().exchange();

            nAtom  = atomStorage.nAtom();
            nGhost = atomStorage.nGhost();
   
            // Check that all atoms are within the processor domain.
            atomStorage.begin(atomIter);
            for ( ; atomIter.notEnd(); ++atomIter) {
               TEST_ASSERT(domain.isInDomain(atomIter->position()));
            }
   
            // Check that all ghosts are outside the processor domain.
            atomStorage.begin(ghostIter);
            for ( ; ghostIter.notEnd(); ++ghostIter) {
               TEST_ASSERT(!domain.isInDomain(ghostIter->position()));
            }

            if (!UTIL_ORTHOGONAL) {
               atomStorage.transformGenToCart(boundary);
            }
            atomStorage.makeSnapshot();
            ++nExchange;
            j = 0;

         } else {

            object().update();

            TEST_ASSERT(nGhost == atomStorage.nGhost());
            TEST_ASSERT(nAtom == atomStorage.nAtom());

            ++ nUpdate;
 
            #if 0
            if (domain.gridRank() == 0) {
               std::cout << "step i = " << i << "  U" << std::endl;
            }
            #endif

         }
   
         TEST_ASSERT(atomStorage.isValid());
         TEST_ASSERT(bondStorage.isValid(atomStorage, domain.communicator(),
                                         true)); 
         #ifdef INTER_ANGLE
         TEST_ASSERT(angleStorage.isValid(atomStorage, 
                     domain.communicator(), true));
         #endif
         #ifdef INTER_DIHEDRAL
         TEST_ASSERT(dihedralStorage.isValid(atomStorage, 
                     domain.communicator(), true));
         #endif

      }
      if (domain.gridRank() == 0) {
         std::cout << std::endl;
         std::cout << "nExchange = " << nExchange << std::endl;
         std::cout << "nUpdate   = " << nUpdate << std::endl;
      }

   }

};

TEST_BEGIN(ExchangerTest)
TEST_ADD(ExchangerTest, testDistribute)
TEST_ADD(ExchangerTest, testAtomExchange)
TEST_ADD(ExchangerTest, testExchange)
TEST_ADD(ExchangerTest, testGhostUpdate)
TEST_ADD(ExchangerTest, testGhostUpdateCycle)
TEST_ADD(ExchangerTest, testExchangeUpdateCycle)
TEST_END(ExchangerTest)

#endif /* EXCHANGER_TEST_H */
