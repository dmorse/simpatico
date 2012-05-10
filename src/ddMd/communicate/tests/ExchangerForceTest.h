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
#include <util/containers/DArray.h>
#include <util/space/Vector.h>

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

   DArray<Vector> forces;

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
      pairPotential.setParamCommunicator(communicator());
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
      pairPotential.readParam(file());
      closeFile();

      object().setPairCutoff(pairPotential.cutoff());
      object().allocate();

      forces.allocate(atomStorage.totalAtomCapacity());

      MaskPolicy policy = MaskBonded;
      std::ifstream configFile("in/config");
      configIo.readConfig(configFile, policy);

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

   void zeroForces()
   {
      AtomIterator  atomIter;
      GhostIterator ghostIter;

      // Zero atom forces
      atomStorage.begin(atomIter);
      for ( ; atomIter.notEnd(); ++atomIter) {
         atomIter->force().zero();
      }

      // Zero ghost forces
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

   void saveForces()
   {
      AtomIterator  atomIter;
      GhostIterator ghostIter;
      int id;

      for (id = 0; id < forces.capacity(); ++id) {
         forces[id].zero();
      }
 
      atomStorage.begin(atomIter);
      for ( ; atomIter.notEnd(); ++atomIter) {
         id = atomIter->id();
         forces[id] = atomIter->force();
      }

      #if 0
      atomStorage.begin(ghostIter);
      for ( ; ghostIter.notEnd(); ++ghostIter) {
         id = ghostIter->id();
         forces[id] = ghostIter->force();
      }
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

   void testInitialForces()
   {
      printMethod(TEST_FUNC);

      int  nAtom  = 0;    // Number of atoms on this processor.
      int  nGhost = 0;    // Number of ghosts on this processor.
      int  nAtomAll  = 0; // Number received on all processors.
      int  myRank = domain.gridRank();

      AtomIterator   atomIter;
      GhostIterator  ghostIter;

      //double range = 0.4;
      //displaceAtoms(range);
      object().exchange();

      // Record number of atoms and ghosts after exchange
      nAtom = atomStorage.nAtom();
      nGhost = atomStorage.nGhost();

      #if 0
      // Update ghost positions
      object().update();

      // Check number of atoms and ghosts on each processor is unchanged.
      TEST_ASSERT(nAtom == atomStorage.nAtom());
      TEST_ASSERT(nGhost == atomStorage.nGhost());
      #endif

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

      // Check that reverse force communication is off (by default)
      // TEST_ASSERT(!pairPotential.forceCommFlag());
      pairPotential.findNeighbors();

      zeroForces();
      pairPotential.setMethodId(2); // N^2 loop
      pairPotential.addForces();
      saveForces();

      zeroForces();
      pairPotential.setMethodId(0); // PairList
      pairPotential.addForces();

      int id;
      std::cout << std::endl;
      atomStorage.begin(atomIter);
      for ( ; atomIter.notEnd(); ++atomIter) {
         id = atomIter->id();
         //std::cout << id << "  "
         //          << forces[id] << "  "
         //          << atomIter->force() << std::endl;
         TEST_ASSERT(eq(forces[id][0], atomIter->force()[0]));
         TEST_ASSERT(eq(forces[id][1], atomIter->force()[1]));
         TEST_ASSERT(eq(forces[id][2], atomIter->force()[2]));
      }

   }

};

TEST_BEGIN(ExchangerForceTest)
TEST_ADD(ExchangerForceTest, testGhostUpdate)
TEST_ADD(ExchangerForceTest, testGhostUpdateCycle)
TEST_ADD(ExchangerForceTest, testInitialForces)
TEST_END(ExchangerForceTest)

#endif /* EXCHANGER_TEST_H */
