#ifndef DDMD_EXCHANGER_FORCE_TEST_H
#define DDMD_EXCHANGER_FORCE_TEST_H

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
#include <ddMd/chemistry/MaskPolicy.h>

#include <util/random/Random.h>
#include <util/containers/DArray.h>
#include <util/space/Vector.h>


#include <ddMd/potentials/pair/PairPotentialImpl.h>
#include <inter/pair/DpdPair.h>

#define TEST_EXCHANGER_FORCE_BOND
#ifdef TEST_EXCHANGER_FORCE_BOND
#include <ddMd/potentials/bond/BondPotentialImpl.h>
#include <inter/bond/HarmonicL0Bond.h>
#ifdef INTER_ANGLE
#include <ddMd/potentials/angle/AnglePotentialImpl.h>
#include <inter/angle/HarmonicAngle.h>
#endif
#ifdef INTER_DIHEDRAL
#include <ddMd/potentials/dihedral/DihedralPotentialImpl.h>
#include <inter/dihedral/MultiHarmonicDihedral.h>
#include <inter/dihedral/CosineDihedral.h>
#endif
#endif 

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

class ExchangerForceTest: public ParamFileTest
{

private:

   Boundary boundary;
   Domain domain;
   Buffer buffer;
   Exchanger exchanger;
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
   int  atomCount;
   bool reverseUpdateFlag;

   DArray<Vector> forces;

   PairPotentialImpl<DpdPair>        pairPotential;
   #ifdef TEST_EXCHANGER_FORCE_BOND
   BondPotentialImpl<HarmonicL0Bond> bondPotential;
   bool hasBondPotential;
   #ifdef INTER_ANGLE
   AnglePotentialImpl<HarmonicAngle> anglePotential;
   bool hasAnglePotential;
   #endif
   #ifdef INTER_DIHEDRAL
   DihedralPotentialImpl<CosineDihedral> dihedralPotential;
   bool hasDihedralPotential;
   #endif
   #endif

public:

   void setUp()
   {
      hasBondPotential = 1;
      #ifdef INTER_ANGLE
      hasAnglePotential = 1;
      #endif
      #ifdef INTER_DIHEDRAL
      hasDihedralPotential = 1;
      #endif
   }

   void initialize()
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
      exchanger.associate(domain, boundary, atomStorage, bondStorage, 
                         #ifdef INTER_ANGLE
                         angleStorage,
                         #endif
                         #ifdef INTER_DIHEDRAL
                         dihedralStorage,
                         #endif
                         buffer);

      pairPotential.setNAtomType(1);
      pairPotential.associate(domain, boundary, atomStorage);
      pairPotential.setReverseUpdateFlag(reverseUpdateFlag);

      #ifdef TEST_EXCHANGER_FORCE_BOND
      if (hasBondPotential) {
         bondPotential.setNBondType(1);
         bondPotential.associate(boundary, bondStorage);
      }
      #ifdef INTER_ANGLE
      if (hasAnglePotential) {
         anglePotential.setNAngleType(1);
         anglePotential.associate(boundary, angleStorage);
      }
      #endif
      #ifdef INTER_DIHEDRAL
      if (hasDihedralPotential) {
         dihedralPotential.setNDihedralType(1);
         dihedralPotential.associate(boundary, dihedralStorage);
      }
      #endif
      #endif

      #ifdef UTIL_MPI
      // Set communicators
      domain.setGridCommunicator(communicator());
      domain.setIoCommunicator(communicator());
      buffer.setIoCommunicator(communicator());
      configIo.setIoCommunicator(communicator());
      random.setIoCommunicator(communicator());
      atomStorage.setIoCommunicator(communicator());
      pairPotential.setIoCommunicator(communicator());
      bondStorage.setIoCommunicator(communicator());
      #ifdef INTER_ANGLE
      angleStorage.setIoCommunicator(communicator());
      #endif
      #ifdef INTER_DIHEDRAL
      dihedralStorage.setIoCommunicator(communicator());
      #endif
      #ifdef TEST_EXCHANGER_FORCE_BOND
      if (hasBondPotential) {
         bondPotential.setIoCommunicator(communicator());
      }
      #ifdef INTER_ANGLE
      if (hasAnglePotential) {
         anglePotential.setIoCommunicator(communicator());
      }
      #endif
      #ifdef INTER_DIHEDRAL
      if (hasDihedralPotential) {
         dihedralPotential.setIoCommunicator(communicator());
      }
      #endif
      #endif // TEST_EXCHANGER_FORCE_BOND
      #else // ifdef UTIL_MPI
      domain.setRank(0);
      #endif

      // Open parameter file
      #ifdef INTER_ANGLE
      #ifdef INTER_DIHEDRAL
      openFile("in/Exchanger_a_d");
      #else // INTER_DIHEDRAL
      openFile("in/Exchanger_a"); 
      #endif // INTER_DIHEDRAL
      #else  // INTER_ANGLE
      openFile("in/Exchanger");
      #endif // INTER_ANGLE

      // Read parameter file
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
      #ifdef TEST_EXCHANGER_FORCE_BOND
      if (hasBondPotential) {
         bondPotential.readParam(file());
      }
      #ifdef INTER_ANGLE
      if (hasAnglePotential) {
         anglePotential.readParam(file());
      }
      #endif
      #ifdef INTER_DIHEDRAL
      if (hasDihedralPotential) {
         dihedralPotential.readParam(file());
      }
      #endif
      #endif
      closeFile();

      exchanger.setPairCutoff(pairPotential.cutoff());
      exchanger.allocate();
      forces.allocate(atomStorage.totalAtomCapacity());

      MaskPolicy policy = MaskBonded;
      //std::ifstream configFile("in/config");
      std::ifstream configFile;
      openInputFile("in/config", configFile);
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

   void displaceAtoms(double range)
   {

      // Set displacement ranges in appropriate coordinate system
      Vector ranges;
      if (UTIL_ORTHOGONAL || atomStorage.isCartesian()) {
         for (int i = 0; i < Dimension; ++i) {
            ranges[i] = range;
         }
      } else {
         for (int i = 0; i < Dimension; ++i) {
            ranges[i] = range/boundary.length(i);
         }
      }
  
      // Displace local atoms
      AtomIterator atomIter;
      double min, max;
      for(int i = 0; i < Dimension; ++i) {
         max = ranges[i];
         min = -max;
         atomStorage.begin(atomIter);
         for ( ; atomIter.notEnd(); ++atomIter) {
            atomIter->position()[i] += random.uniform(min, max);
         }
      }

   }

   void zeroForces()
   {
      // Zero atom forces
      AtomIterator  atomIter;
      atomStorage.begin(atomIter);
      for ( ; atomIter.notEnd(); ++atomIter) {
         atomIter->force().zero();
      }

      // Zero ghost forces
      if (atomStorage.nGhost()) {
         GhostIterator ghostIter;
         atomStorage.begin(ghostIter);
         for ( ; ghostIter.notEnd(); ++ghostIter) {
            ghostIter->force().zero();
         }
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

   void computeForces() {
      zeroForces();
      pairPotential.computeForces();
      #ifdef TEST_EXCHANGER_FORCE_BOND
      if (hasBondPotential) {
         bondPotential.computeForces();
      }
      #ifdef INTER_ANGLE
      if (hasAnglePotential) {
         anglePotential.computeForces();
      }
      #endif
      #ifdef INTER_DIHEDRAL
      if (hasDihedralPotential) {
         dihedralPotential.computeForces();
      }
      #endif
      #endif
    }
   void saveForces()
   {
      int id;
      for (id = 0; id < forces.capacity(); ++id) {
         forces[id].zero();
      }
 
      AtomIterator atomIter;
      atomStorage.begin(atomIter);
      for ( ; atomIter.notEnd(); ++atomIter) {
         id = atomIter->id();
         forces[id] = atomIter->force();
      }

      #if 0
      GhostIterator ghostIter;
      atomStorage.begin(ghostIter);
      for ( ; ghostIter.notEnd(); ++ghostIter) {
         id = ghostIter->id();
         forces[id] = ghostIter->force();
      }
      #endif
   }

   void exchangeNotify() {
      bondStorage.unsetNTotal();
      #ifdef INTER_ANGLE
      angleStorage.unsetNTotal();
      #endif
      #ifdef INTER_DIHEDRAL
      dihedralStorage.unsetNTotal();
      #endif
   }


   void testGhostUpdateF() {
      printMethod(TEST_FUNC);
      reverseUpdateFlag = false;
      testGhostUpdate();
   }

   void testGhostUpdateR() {
      printMethod(TEST_FUNC);
      reverseUpdateFlag = true;
      testGhostUpdate();
   }

   void testGhostUpdate()
   {
      initialize();

      int  nAtom  = 0;    // Number of atoms on this processor.
      int  nGhost = 0;    // Number of ghosts on this processor.
      int  nAtomAll  = 0; // Number received on all processors.
      int  myRank = domain.gridRank();

      double range = 0.4;
      displaceAtoms(range);
      exchanger.exchange();
      exchangeNotify();

      // Record number of atoms and ghosts after exchange
      nAtom = atomStorage.nAtom();
      nGhost = atomStorage.nGhost();

      // Transform to Cartesian coordinates
      if (!UTIL_ORTHOGONAL) {
         atomStorage.transformGenToCart(boundary);
      }

      // Update ghost positions
      exchanger.update();

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

      // Transform to generalized coordinates
      if (!UTIL_ORTHOGONAL) {
         atomStorage.transformCartToGen(boundary);
      }

      // Check that all atoms are within the processor domain.
      AtomIterator   atomIter;
      atomStorage.begin(atomIter);
      for ( ; atomIter.notEnd(); ++atomIter) {
         TEST_ASSERT(domain.isInDomain(atomIter->position()));
      }

      // Check that all ghosts are outside the processor domain.
      GhostIterator  ghostIter;
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

   void testGhostUpdateCycleF() {
      printMethod(TEST_FUNC);
      reverseUpdateFlag = false;
      testGhostUpdateCycle();
   }

   void testGhostUpdateCycleR() {
      printMethod(TEST_FUNC);
      reverseUpdateFlag = true;
      testGhostUpdateCycle();
   }

   void testGhostUpdateCycle()
   {
      initialize();

      int  nAtom  = 0;    // Number of atoms on this processor.
      int  nGhost = 0;    // Number of ghosts on this processor.

      AtomIterator   atomIter;
      GhostIterator  ghostIter;

      double range = 0.05;
      //displaceAtoms(range);

      exchanger.exchange();
      exchangeNotify();
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

      range = 0.05;
      for (int i=0; i < 4; ++i) {

         TEST_ASSERT(atomStorage.isCartesian());
         displaceAtoms(range);

         for (int j=0; j < 2; ++j) {
            exchanger.update();
            TEST_ASSERT(nGhost == atomStorage.nGhost());
            TEST_ASSERT(nAtom == atomStorage.nAtom());
            displaceAtoms(range);
         }

         // Transform to generalized coordinates
         if (!UTIL_ORTHOGONAL) {
            atomStorage.transformCartToGen(boundary);
         }

         atomStorage.clearSnapshot();
         exchanger.exchange();
         exchangeNotify();

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

   void testInitialForcesF() {
      printMethod(TEST_FUNC);
      reverseUpdateFlag = false;
      testInitialForces();
   }

   void testInitialForcesR() {
      printMethod(TEST_FUNC);
      reverseUpdateFlag = true;
      testInitialForces();
   }

   void testInitialForces()
   {
      initialize();

      int  nAtom  = 0;    // Number of atoms on this processor.
      int  nGhost = 0;    // Number of ghosts on this processor.

      atomStorage.clearSnapshot();
      exchanger.exchange();
      exchangeNotify();

      // Record number of atoms and ghosts after exchange
      nAtom = atomStorage.nAtom();
      nGhost = atomStorage.nGhost();

      // Check that all atoms are within the processor domain.
      AtomIterator  atomIter;
      atomStorage.begin(atomIter);
      for ( ; atomIter.notEnd(); ++atomIter) {
         TEST_ASSERT(domain.isInDomain(atomIter->position()));
      }

      // Check that all ghosts are outside the processor domain.
      GhostIterator ghostIter;
      atomStorage.begin(ghostIter);
      for ( ; ghostIter.notEnd(); ++ghostIter) {
         TEST_ASSERT(!domain.isInDomain(ghostIter->position()));
      }

      // Build Cell and Pair lists
      pairPotential.buildCellList();
      if (!UTIL_ORTHOGONAL) {
         atomStorage.transformGenToCart(boundary);
      }
      pairPotential.buildPairList();
      atomStorage.makeSnapshot();

      // Update ghost positions
      exchanger.update();

      // Check number of atoms and ghosts on each processor is unchanged.
      TEST_ASSERT(nAtom == atomStorage.nAtom());
      TEST_ASSERT(nGhost == atomStorage.nGhost());

      // Check that all atoms are accounted for after atom and ghost exchanges.
      nAtom = atomStorage.nAtom();
      int  nAtomAll  = 0; // Number received on all processors.
      communicator().Reduce(&nAtom, &nAtomAll, 1, MPI::INT, MPI::SUM, 0);
      int  myRank = domain.gridRank();
      if (myRank == 0) {
         // std::cout << "Total atom count (post ghost exchange) = " 
         //           << nAtomAll << std::endl;
         TEST_ASSERT(nAtomAll == atomCount);
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

      TEST_ASSERT(pairPotential.reverseUpdateFlag() == reverseUpdateFlag);

      // Compute forces etc. with N^2 loop (MethodId = 2)
      pairPotential.setMethodId(2);
      computeForces();
      if (reverseUpdateFlag) {
         exchanger.reverseUpdate();
      }
      saveForces();
      int nPairNSq;
      pairPotential.computeNPair(domain.communicator());
      if (domain.communicator().Get_rank() == 0) {
         nPairNSq = pairPotential.nPair();
      }
      double pairEnergyNSq;
      pairPotential.unsetEnergy();
      pairPotential.computeEnergy(domain.communicator());
      if (domain.communicator().Get_rank() == 0) {
         pairEnergyNSq = pairPotential.energy();
      }

      // Compute forces etc. with pair list (MethodId = 0)
      pairPotential.setMethodId(0); 
      computeForces();
      if (reverseUpdateFlag) {
         exchanger.reverseUpdate();
      }
      int nPairList;
      pairPotential.computeNPair(domain.communicator());
      if (domain.communicator().Get_rank() == 0) {
         nPairList = pairPotential.nPair();
      }
      double energyList;
      pairPotential.unsetEnergy();
      pairPotential.computeEnergy(domain.communicator());
      if (domain.communicator().Get_rank() == 0) {
         energyList = pairPotential.energy();
      }
      Tensor pairStress;
      pairPotential.unsetStress();
      pairPotential.computeStress(domain.communicator());
      if (domain.communicator().Get_rank() == 0) {
         pairStress = pairPotential.stress();
      }
      #ifdef TEST_EXCHANGER_FORCE_BOND
      Tensor bondStress;
      if (hasBondPotential) {
         bondPotential.computeStress(domain.communicator());
         if (domain.communicator().Get_rank() == 0) {
            bondStress = bondPotential.stress();
         }
      }
      #ifdef INTER_ANGLE
      Tensor angleStress;
      if (hasDihedralPotential) {
         anglePotential.unsetStress();
         anglePotential.computeStress(domain.communicator());
         if (domain.communicator().Get_rank() == 0) {
            angleStress = anglePotential.stress();
         }
      }
      #endif
      #ifdef INTER_ANGLE
      Tensor dihedralStress;
      if (hasDihedralPotential) {
         dihedralPotential.unsetStress();
         dihedralPotential.computeStress(domain.communicator());
         if (domain.communicator().Get_rank() == 0) {
            dihedralStress = dihedralPotential.stress();
         }
      }
      #endif
      #endif

      // Compare Nsq and PairList values of nPair and Energy
      if (domain.communicator().Get_rank() == 0) {
         TEST_ASSERT(nPairNSq == nPairList);
         TEST_ASSERT(eq(pairEnergyNSq, energyList));
      }

      // Compare NSq and pairlist forces, accumulate total
      Vector totForce;
      Vector nodeForce;
      int id, i;
      bool isEqual; 
      nodeForce.zero();
      atomStorage.begin(atomIter);
      for ( ; atomIter.notEnd(); ++atomIter) {
         id = atomIter->id();
         for (int i = 0; i < 3; ++i) {
            isEqual = eq(forces[id][i], atomIter->force()[i]);
            if (!isEqual) {
               std::cout << id << "  "
                         << forces[id][i] << "    "
                         << atomIter->force()[i] << std::endl;

            }
            TEST_ASSERT(isEqual);
         }
         nodeForce += atomIter->force();
      }

      // Check that total force is zero (on master node)
      communicator().Reduce(&nodeForce[0], &totForce[0], 3, MPI::DOUBLE, MPI::SUM, 0);
      if (communicator().Get_rank() == 0) {
         TEST_ASSERT(eq(totForce[0], 0.0));
         TEST_ASSERT(eq(totForce[1], 0.0));
         TEST_ASSERT(eq(totForce[2], 0.0));
      }

      // Test computeForcesAndStress methods 
      pairPotential.setMethodId(0); // PairList
      zeroForces();
      pairPotential.unsetStress();
      pairPotential.computeForcesAndStress(domain.communicator());
      #ifdef TEST_EXCHANGER_FORCE_BOND
      if (hasBondPotential) {
         bondPotential.unsetStress();
         bondPotential.computeForcesAndStress(domain.communicator());
      }
      #ifdef INTER_ANGLE
      if (hasAnglePotential) {
         anglePotential.unsetStress();
         anglePotential.computeForcesAndStress(domain.communicator());
      }
      #endif
      #ifdef INTER_DIHEDRAL
      if (hasDihedralPotential) {
         dihedralPotential.unsetStress();
         dihedralPotential.computeForcesAndStress(domain.communicator());
      }
      #endif
      #endif TEST_EXCHANGER_FORCE_BOND
      if (reverseUpdateFlag) {
         exchanger.reverseUpdate();
      }

      if (domain.communicator().Get_rank() == 0) {
         Tensor pairStress2 = pairPotential.stress();
         #ifdef TEST_EXCHANGER_FORCE_BOND
         Tensor bondStress2;
         if (hasBondPotential) {
            bondStress2 = bondPotential.stress();
         }
         #endif
         for (int i = 0; i < Dimension; ++i) {
            for (int j = 0; j < Dimension; ++j) {
               TEST_ASSERT(eq(pairStress(i, j), pairStress2(i, j)));
               #ifdef TEST_EXCHANGER_FORCE_BOND
               if (hasBondPotential) {
                  TEST_ASSERT(eq(bondStress(i, j), bondStress2(i, j)));
               }
               #endif // TEST_EXCHANGER_FORCE_BOND
            }
         }
      }
      atomStorage.begin(atomIter);
      for ( ; atomIter.notEnd(); ++atomIter) {
         id = atomIter->id();
         TEST_ASSERT(eq(forces[id][0], atomIter->force()[0]));
         TEST_ASSERT(eq(forces[id][1], atomIter->force()[1]));
         TEST_ASSERT(eq(forces[id][2], atomIter->force()[2]));
      }

   }

   void testForceCycleF() {
      printMethod(TEST_FUNC);
      reverseUpdateFlag = true;
      testForceCycle();
   }

   void testForceCycleR() {
      printMethod(TEST_FUNC);
      reverseUpdateFlag = false;
      testForceCycle();
   }

   void testForceCycle()
   {
      initialize();

      int  nAtom  = 0;    // Number of atoms on this processor.
      int  nGhost = 0;    // Number of ghosts on this processor.
      bool needExchange;

      TEST_ASSERT(pairPotential.reverseUpdateFlag() == reverseUpdateFlag);

      // double range = 0.1;
      // displaceAtoms(range);

      atomStorage.clearSnapshot();
      if (!UTIL_ORTHOGONAL) {
         TEST_ASSERT(!atomStorage.isCartesian());
      }
      exchanger.exchange();
      exchangeNotify();

      nAtom = atomStorage.nAtom();
      nGhost = atomStorage.nGhost();

      // Assert that all atoms are within the processor domain.
      AtomIterator   atomIter;
      atomStorage.begin(atomIter);
      for ( ; atomIter.notEnd(); ++atomIter) {
         TEST_ASSERT(domain.isInDomain(atomIter->position()));
      }

      // Assert that all ghosts are outside the processor domain.
      GhostIterator  ghostIter;
      atomStorage.begin(ghostIter);
      for ( ; ghostIter.notEnd(); ++ghostIter) {
         TEST_ASSERT(!domain.isInDomain(ghostIter->position()));
      }

      TEST_ASSERT(atomStorage.isValid());
      TEST_ASSERT(bondStorage.isValid(atomStorage, domain.communicator(), 
                  true));

      // Calculate forces with PairList
      TEST_ASSERT(!atomStorage.isCartesian());
      //pairPotential.findNeighbors();
      pairPotential.buildCellList();
      if (!UTIL_ORTHOGONAL) {
         atomStorage.transformGenToCart(boundary);
      }
      TEST_ASSERT(atomStorage.isCartesian());
      pairPotential.buildPairList();
      atomStorage.makeSnapshot();
      TEST_ASSERT(atomStorage.isCartesian());

      // Compute forces (Pairlist method)
      pairPotential.setMethodId(0); 
      computeForces();
      if (reverseUpdateFlag) {
         exchanger.reverseUpdate();
      }

      double pairEnergyNSq, energyList, energyF;
      int    nPairNSq, nPairList, nPairF;

      double range = 0.02;

      int i = 0;  // step counter
      int j = 0;  // exchange counter
      for ( ; i < 100; ++i) {

         TEST_ASSERT(atomStorage.isCartesian());
         displaceAtoms(range);
   
         // Check if exchange and reneighboring is necessary
         needExchange = atomStorage.needExchange(domain.communicator(), 
                                                 pairPotential.skin());
         if (needExchange) {
            //if (domain.isMaster()) {
            //   std::cout << "Step " << i << ",  exchange " << j << std::endl;
            //}
            atomStorage.clearSnapshot();
            TEST_ASSERT(atomStorage.isCartesian());
            if (!UTIL_ORTHOGONAL) {
               atomStorage.transformCartToGen(boundary);
               TEST_ASSERT(!atomStorage.isCartesian());
            }

            exchanger.exchange();
            exchangeNotify();

            // Confirm that all atoms are within the processor domain.
            atomStorage.begin(atomIter);
            for ( ; atomIter.notEnd(); ++atomIter) {
               TEST_ASSERT(domain.isInDomain(atomIter->position()));
            }
   
            // Confirm that all ghosts are outside the processor domain.
            atomStorage.begin(ghostIter);
            for ( ; ghostIter.notEnd(); ++ghostIter) {
               TEST_ASSERT(!domain.isInDomain(ghostIter->position()));
            }

            pairPotential.buildCellList();
            if (!UTIL_ORTHOGONAL) {
               TEST_ASSERT(!atomStorage.isCartesian());
               atomStorage.transformGenToCart(boundary);
            } else {
               TEST_ASSERT(atomStorage.isCartesian());
            }
            pairPotential.buildPairList();
            atomStorage.makeSnapshot();
            TEST_ASSERT(atomStorage.isCartesian());

            ++j;

         } else {

            exchanger.update();
         }
   
         nAtom  = atomStorage.nAtom();
         nGhost = atomStorage.nGhost();

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

         // Calculate forces by an N^2 loop.
         computeForces();
         if (reverseUpdateFlag) {
            exchanger.reverseUpdate();
         }
         saveForces();
         pairPotential.computeNPair(domain.communicator());
         if (domain.communicator().Get_rank() == 0) {
            nPairNSq = pairPotential.nPair();
         }
         pairPotential.computeEnergy(domain.communicator());
         if (domain.communicator().Get_rank() == 0) {
            pairEnergyNSq = pairPotential.energy();
         }

         // Calculate forces etc. via pair list.   
         pairPotential.setMethodId(0); // PairList
         computeForces();
         if (reverseUpdateFlag) {
            exchanger.reverseUpdate();
         }
         pairPotential.computeEnergy(domain.communicator());
         if (domain.communicator().Get_rank() == 0) {
            energyList = pairPotential.energy();
         }
         pairPotential.computeNPair(domain.communicator());
         if (domain.communicator().Get_rank() == 0) {
            nPairList = pairPotential.nPair();
         }
      
         Vector nodeForce;
         nodeForce.zero();
         int id;
         atomStorage.begin(atomIter);
         for ( ; atomIter.notEnd(); ++atomIter) {
            id = atomIter->id();
            TEST_ASSERT(eq(forces[id][0], atomIter->force()[0]));
            TEST_ASSERT(eq(forces[id][1], atomIter->force()[1]));
            TEST_ASSERT(eq(forces[id][2], atomIter->force()[2]));
            nodeForce += atomIter->force();
         }
   
         // Check that total force is zero, different methods agree.
         Vector totForce;
         communicator().Reduce(&nodeForce[0], &totForce[0], 3, MPI::DOUBLE, MPI::SUM, 0);
         if (communicator().Get_rank() == 0) {
            TEST_ASSERT(eq(totForce[0], 0.0));
            TEST_ASSERT(eq(totForce[1], 0.0));
            TEST_ASSERT(eq(totForce[2], 0.0));
            TEST_ASSERT(nPairNSq == nPairList);
            TEST_ASSERT(eq(pairEnergyNSq, energyList));
         }

         if (reverseUpdateFlag && needExchange) {

            // Calculate forces via pair list, without reverse communication
            zeroForces();
            pairPotential.setReverseUpdateFlag(false); 
            TEST_ASSERT(atomStorage.isCartesian());
            if (!UTIL_ORTHOGONAL) {
               atomStorage.transformCartToGen(boundary);
            }
            TEST_ASSERT(!atomStorage.isCartesian());
            pairPotential.buildCellList();
            if (!UTIL_ORTHOGONAL) {
               atomStorage.transformGenToCart(boundary);
            }
            pairPotential.buildPairList();
            TEST_ASSERT(atomStorage.isCartesian());
            pairPotential.setMethodId(0);    
            computeForces();
            saveForces();
            pairPotential.computeEnergy(domain.communicator());
            if (domain.communicator().Get_rank() == 0) {
               energyF = pairPotential.energy();
            }
            pairPotential.computeNPair(domain.communicator());
            if (domain.communicator().Get_rank() == 0) {
               nPairF = pairPotential.nPair();
            }
   
            // Compare values for atom forces
            atomStorage.begin(atomIter);
            for ( ; atomIter.notEnd(); ++atomIter) {
               id = atomIter->id();
               TEST_ASSERT(eq(forces[id][0], atomIter->force()[0]));
               TEST_ASSERT(eq(forces[id][1], atomIter->force()[1]));
               TEST_ASSERT(eq(forces[id][2], atomIter->force()[2]));
            }

            // Compare energy and nPair
            if (communicator().Get_rank() == 0) {
               TEST_ASSERT(nPairF == nPairList);
               TEST_ASSERT(eq(energyF, energyList));
            }

            // Reset: Recompute neighbor list for use with reverse communication
            pairPotential.setReverseUpdateFlag(true); 
            TEST_ASSERT(atomStorage.isCartesian());
            atomStorage.transformCartToGen(boundary);
            TEST_ASSERT(!atomStorage.isCartesian());
            //pairPotential.findNeighbors(); 
            pairPotential.buildCellList();
            if (!UTIL_ORTHOGONAL) {
               atomStorage.transformGenToCart(boundary);
            }
            pairPotential.buildPairList();

         }

      }

      #if 0
      if (domain.isMaster()) {
         std::cout << std::endl;
         std::cout << "Step " << i << ",  exchange " << j << std::endl;
      }
      #endif

   }

};

TEST_BEGIN(ExchangerForceTest)
TEST_ADD(ExchangerForceTest, testGhostUpdateF)
TEST_ADD(ExchangerForceTest, testGhostUpdateR)
TEST_ADD(ExchangerForceTest, testGhostUpdateCycleF)
TEST_ADD(ExchangerForceTest, testGhostUpdateCycleR)
TEST_ADD(ExchangerForceTest, testInitialForcesF)
TEST_ADD(ExchangerForceTest, testInitialForcesR)
TEST_ADD(ExchangerForceTest, testForceCycleF)
TEST_ADD(ExchangerForceTest, testForceCycleR)
TEST_END(ExchangerForceTest)

#endif /* EXCHANGER_TEST_H */
