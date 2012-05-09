#ifndef DDMD_EXCHANGER_FORCE_TEST_H
#define DDMD_EXCHANGER_FORCE_TEST_H

#include <ddMd/configIos/ConfigIo.h>
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
#include <ddMd/chemistry/MaskPolicy.h>
#include <ddMd/potentials/pair/PairPotentialImpl.h>

#include <inter/pair/DpdPair.h>

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
using namespace Inter;
using namespace DdMd;

class ExchangerForceTest: public ParamFileTest<Exchanger>
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
   ConfigIo configIo;
   Random random;

   PairPotentialImpl<DpdPair> pairPotential;

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

      pairPotential.setNAtomType(1);
      pairPotential.associate(domain, boundary, atomStorage);

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

      // Read parameter file
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
      //pairPotential.readParam(file());
      closeFile();

      MaskPolicy policy = MaskBonded;
      std::ifstream configFile("in/config");
      configIo.readConfig(configFile, policy);

      object().allocate();
      object().setPairCutoff(0.5);

      int  nAtom = 0;     // Number received on this processor.
      int  nAtomAll  = 0; // Number received on all processors.

      // Check that all atoms are accounted for after distribution.
      nAtom = atomStorage.nAtom();
      communicator().Reduce(&nAtom, &nAtomAll, 1, MPI::INT, MPI::SUM, 0);
      if (domain.gridRank() == 0) {
         //std::cout << std::endl;
         // std::cout << "Total atom count (post-distribute) = " 
         //          << nAtomAll << std::endl;
         atomCount = nAtomAll;
      }

   }

   void displaceAtoms(double max)
   {
      double min = -max;
      AtomIterator atomIter;
      for(int i = 0; i < 3; i++) {
         atomStorage.begin(atomIter);
         for ( ; atomIter.notEnd(); ++atomIter) {
            atomIter->position()[i] += random.uniform(min, max);
         }
      }
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

      double range = 0.4;
      displaceAtoms(range);
      object().exchange();

      // Record number of atoms and ghosts after exchange
      nAtom = atomStorage.nAtom();
      nGhost = atomStorage.nGhost();

      // Update ghost positions
      object().update();

      // Check number of atoms and ghosts on each processor is unchanged.
      TEST_ASSERT(nAtom == atomStorage.nAtom());
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
      int  nAtomAll  = 0; // Number received on all processors.
      int  myRank = domain.gridRank();

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

      range = 0.1;
      for (int i=0; i < 3; ++i) {

         displaceAtoms(range);

         for (int j=0; j < 3; ++j) {
            object().update();
            TEST_ASSERT(nGhost == atomStorage.nGhost());
            TEST_ASSERT(nAtom == atomStorage.nAtom());
            displaceAtoms(range);
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
      }

   }


   void zeroForces()
   {
      AtomIterator  atomIter;
      GhostIterator ghostIter;

      atomStorage.begin(atomIter);
      for ( ; atomIter.notEnd(); ++atomIter) {
         atomIter->force().zero();
      }
      atomStorage.begin(ghostIter);
      for ( ; ghostIter.notEnd(); ++ghostIter) {
         ghostIter->force().zero();
      }

   }

   void writeForces()
   {
      std::cout << std::endl;

      AtomIterator  atomIter;
      atomStorage.begin(atomIter);
      for ( ; atomIter.notEnd(); ++atomIter) {
         std::cout << atomIter->force() << "  " << 0 << std::endl;
      }

      #if 0
      GhostIterator ghostIter;
      atomStorage.begin(ghostIter);
      for ( ; ghostIter.notEnd(); ++ghostIter) {
         std::cout << ghostIter->force() << "  " << 1 << std::endl;
      }
      #endif

   }

   #if 0
   void saveForces()
   {
      AtomIterator  atomIter;
      GhostIterator ghostIter;

      atomStorage.begin(atomIter);
      for ( ; atomIter.notEnd(); ++atomIter) {
         force_[i] = atomIter->force();
      }
      atomStorage.begin(ghostIter);
      for ( ; ghostIter.notEnd(); ++ghostIter) {
         force_[i] = ghostIter->force();
      }

   }
   #endif

   void testRead1()
   {
      printMethod(TEST_FUNC);
      //readAtoms("in/positions1");
     
      TEST_ASSERT(atomStorage.nAtom()  == 3); 
      TEST_ASSERT(atomStorage.nGhost() == 0); 

      AtomIterator iter;
      std::cout << std::endl;
      atomStorage.begin(iter);
      for ( ; iter.notEnd(); ++iter) {
         std::cout << iter->position() << std::endl;
      }      

      pairPotential.setMethodId(2); // N^2 loop
      pairPotential.addForces();

      std::cout << std::endl;
      atomStorage.begin(iter);
      for ( ; iter.notEnd(); ++iter) {
         std::cout << iter->force() << std::endl;
      }      
   }

};

TEST_BEGIN(ExchangerForceTest)
TEST_ADD(ExchangerForceTest, testGhostUpdate)
TEST_ADD(ExchangerForceTest, testGhostUpdateCycle)
TEST_END(ExchangerForceTest)

#endif /* EXCHANGER_TEST_H */
