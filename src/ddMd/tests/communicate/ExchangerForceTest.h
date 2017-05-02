#ifndef DDMD_EXCHANGER_FORCE_TEST_H
#define DDMD_EXCHANGER_FORCE_TEST_H

#include <ddMd/configIos/DdMdConfigIo.h>
#include <ddMd/communicate/Domain.h>
#include <ddMd/communicate/Buffer.h>
#include <ddMd/communicate/Exchanger.h>
#include <ddMd/communicate/GroupDistributor.h>
#include <ddMd/communicate/GroupDistributor.tpp>
#include <ddMd/communicate/GroupCollector.h>
#include <ddMd/communicate/GroupCollector.tpp>
#include <ddMd/storage/AtomStorage.h>
#include <ddMd/storage/AtomIterator.h>
#include <ddMd/storage/GhostIterator.h>
#include <ddMd/storage/BondStorage.h>
#ifdef SIMP_ANGLE
#include <ddMd/storage/AngleStorage.h>
#endif
#ifdef SIMP_DIHEDRAL
#include <ddMd/storage/DihedralStorage.h>
#endif
#include <ddMd/chemistry/MaskPolicy.h>

#include <util/random/Random.h>
#include <util/containers/DArray.h>
#include <util/space/Vector.h>

#include <ddMd/potentials/pair/PairPotentialImpl.h>
#include <simp/interaction/pair/DpdPair.h>

#include <ddMd/potentials/bond/BondPotentialImpl.h>
#include <simp/interaction/bond/HarmonicL0Bond.h>
#ifdef SIMP_ANGLE
#include <ddMd/potentials/angle/AnglePotentialImpl.h>
#include <simp/interaction/angle/HarmonicAngle.h>
#endif
#ifdef SIMP_DIHEDRAL
#include <ddMd/potentials/dihedral/DihedralPotentialImpl.h>
#include <simp/interaction/dihedral/MultiHarmonicDihedral.h>
#include <simp/interaction/dihedral/CosineDihedral.h>
#endif

#ifdef UTIL_MPI
#ifndef TEST_MPI
#define TEST_MPI
#endif
#endif

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>
#include <test/ParamFileTest.h>

using namespace Util;
using namespace Simp;
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
   #ifdef SIMP_ANGLE
   AngleStorage angleStorage;
   #endif
   #ifdef SIMP_DIHEDRAL
   DihedralStorage dihedralStorage;
   #endif
   DdMdConfigIo configIo;
   Random random;
   int  atomCount;
   bool reverseUpdateFlag;
   bool hasAngles;
   bool hasDihedrals;

   DArray<Vector> forces;

   PairPotentialImpl<DpdPair> pairPotential;
   BondPotentialImpl<HarmonicL0Bond> bondPotential;
   #ifdef SIMP_ANGLE
   AnglePotentialImpl<HarmonicAngle> anglePotential;
   #endif
   #ifdef SIMP_DIHEDRAL
   DihedralPotentialImpl<CosineDihedral> dihedralPotential;
   #endif

   void initialize();
   void displaceAtoms(double range);
   void zeroForces();
   void computeForces();
   void writeForces();
   void saveForces();
   bool isExchangeNeeded(double skin);
   void exchangeNotify();

   void testGhostUpdate();
   void testGhostUpdateCycle();
   void testInitialForces();
   void testForceCycle();

public:

   void setUp();
   void testGhostUpdateF();
   void testGhostUpdateR();
   void testGhostUpdateCycleF();
   void testGhostUpdateCycleR();
   void testInitialForcesF();
   void testInitialForcesR();
   void testForceCycleF();
   void testForceCycleR();
};

void ExchangerForceTest::setUp()
{
   hasAngles = 0;
   hasDihedrals = 0;
}

void ExchangerForceTest::initialize()
{

   // Set connections between atomDistributors
   domain.setBoundary(boundary);
   atomStorage.associate(domain, boundary, buffer);
   bondStorage.associate(domain, atomStorage, buffer);
   exchanger.associate(domain, boundary, atomStorage, buffer);
   exchanger.addGroupExchanger(bondStorage);
   #ifdef SIMP_ANGLE
   angleStorage.associate(domain, atomStorage, buffer);
   if (hasAngle) {
      exchanger.addGroupExchanger(angleStorage);
   }
   #endif
   #ifdef SIMP_DIHEDRAL
   dihedralStorage.associate(domain, atomStorage, buffer);
   if (hasDihedral) {
      exchanger.addGroupExchanger(dihedralStorage);
   }
   #endif
   configIo.associate(domain, boundary, atomStorage, bondStorage, 
                      #ifdef SIMP_ANGLE
                      angleStorage,
                      #endif
                      #ifdef SIMP_DIHEDRAL
                      dihedralStorage,
                      #endif
                      buffer);

   #ifdef UTIL_MPI
   // Set communicators
   domain.setGridCommunicator(communicator());
   domain.setIoCommunicator(communicator());
   buffer.setIoCommunicator(communicator());
   configIo.setIoCommunicator(communicator());
   random.setIoCommunicator(communicator());
   atomStorage.setIoCommunicator(communicator());
   bondStorage.setIoCommunicator(communicator());
   #ifdef SIMP_ANGLE
   angleStorage.setIoCommunicator(communicator());
   #endif
   #ifdef SIMP_DIHEDRAL
   dihedralStorage.setIoCommunicator(communicator());
   #endif
   #else // ifdef UTIL_MPI
   domain.setRank(0);
   #endif

   // Open parameter file
   if (hasDihedrals) {
      openFile("in/Exchanger_a_d");
   } else 
   if (hasAngles) {
      openFile("in/Exchanger_a"); 
   } else {
      openFile("in/Exchanger");
   }

   // Read parameter file
   domain.readParam(file());
   buffer.readParam(file());
   random.readParam(file());

   atomStorage.readParam(file());
   #ifdef SIMP_BOND
   bondStorage.readParam(file());
   #endif
   #ifdef SIMP_ANGLE
   if (hasAngle) {
      angleStorage.readParam(file());
   }
   #endif
   #ifdef SIMP_DIHEDRAL
   if (hasDihedral) {
      dihedralStorage.readParam(file());
   }
   #endif
   
   pairPotential.associate(domain, boundary, atomStorage);
   pairPotential.setIoCommunicator(communicator());
   pairPotential.setNAtomType(1);
   pairPotential.setReverseUpdateFlag(reverseUpdateFlag);
   pairPotential.readParam(file());

   bondPotential.associate(boundary, bondStorage);
   bondPotential.setIoCommunicator(communicator());
   bondPotential.setNBondType(1);
   bondPotential.readParam(file());

   #ifdef SIMP_ANGLE
   if (hasAngles) {
      anglePotential.setIoCommunicator(communicator());
      anglePotential.associate(boundary, angleStorage);
      anglePotential.setNAngleType(1);
      anglePotential.readParam(file());
   }
   #endif
   #ifdef SIMP_DIHEDRAL
   if (hasDihedrals) {
      dihedralPotential.associate(boundary, dihedralStorage);
      dihedralPotential.setIoCommunicator(communicator());
      dihedralPotential.setNDihedralType(1);
      dihedralPotential.readParam(file());
   }
   #endif

   closeFile();

   exchanger.setPairCutoff(pairPotential.cutoff());
   exchanger.allocate();
   forces.allocate(atomStorage.totalAtomCapacity());

   MaskPolicy policy = MaskBonded;
   std::ifstream configFile;
   openInputFile("in/config", configFile);
   configIo.readConfig(configFile, policy);

   int  nAtom = 0;    // Number received on this processor.
   int  nAtomAll = 0; // Number received on all processors.

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

void ExchangerForceTest::displaceAtoms(double range)
{

   // Set displacement ranges in appropriate coordinate system
   // Input parameter is range in Cartesian coordinates
   Vector ranges;
   if (atomStorage.isCartesian()) {
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
   for (atomStorage.begin(atomIter); atomIter.notEnd(); ++atomIter) {
      for (int i = 0; i < Dimension; ++i) {
         max = ranges[i];
         min = -max;
         atomIter->position()[i] += random.uniform(min, max);
      }
   }

}

void ExchangerForceTest::zeroForces()
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

void ExchangerForceTest::writeForces()
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

void ExchangerForceTest::computeForces() {
   zeroForces();
   pairPotential.computeForces();
   bondPotential.computeForces();
   #ifdef SIMP_ANGLE
   if (hasAngles) {
      anglePotential.computeForces();
   }
   #endif
   #ifdef SIMP_DIHEDRAL
   if (hasDihedrals) {
      dihedralPotential.computeForces();
   }
   #endif
}

void ExchangerForceTest::saveForces()
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

void ExchangerForceTest::exchangeNotify() {
   bondStorage.unsetNTotal();
   #ifdef SIMP_ANGLE
   if (hasAngles) {
      angleStorage.unsetNTotal();
   }
   #endif
   #ifdef SIMP_DIHEDRAL
   if (hasDihedrals) {
      dihedralStorage.unsetNTotal();
   }
   #endif
}


void ExchangerForceTest::testGhostUpdateF() {
   printMethod(TEST_FUNC);
   reverseUpdateFlag = false;
   testGhostUpdate();
}

void ExchangerForceTest::testGhostUpdateR() {
   printMethod(TEST_FUNC);
   reverseUpdateFlag = true;
   testGhostUpdate();
}

void ExchangerForceTest::testGhostUpdate()
{
   initialize();

   int  nAtom  = 0;    // Number of atoms on this processor.
   int  nGhost = 0;    // Number of ghosts on this processor.
   int  nAtomAll = 0; // Number received on all processors.
   int  myRank = domain.gridRank();

   double range = 0.2;
   displaceAtoms(range);
   atomStorage.clearSnapshot();
   exchanger.exchange();
   exchangeNotify();

   // Record number of atoms and ghosts after exchange
   nAtom = atomStorage.nAtom();
   nGhost = atomStorage.nGhost();

   // Transform to Cartesian coordinates
   atomStorage.transformGenToCart(boundary);
   atomStorage.makeSnapshot();

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
   atomStorage.transformCartToGen(boundary);

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
   #ifdef SIMP_ANGLE
   if (hasAngles) {
      TEST_ASSERT(angleStorage.isValid(atomStorage, 
                  domain.communicator(), true));
   }
   #endif
   #ifdef SIMP_DIHEDRAL
   if (hasDihedrals) {
      TEST_ASSERT(dihedralStorage.isValid(atomStorage, 
                  domain.communicator(), true));
   }
   #endif

}

void ExchangerForceTest::testGhostUpdateCycleF() {
   printMethod(TEST_FUNC);
   reverseUpdateFlag = false;
   testGhostUpdateCycle();
}

void ExchangerForceTest::testGhostUpdateCycleR() {
   printMethod(TEST_FUNC);
   reverseUpdateFlag = true;
   testGhostUpdateCycle();
}

void ExchangerForceTest::testGhostUpdateCycle()
{
   initialize();

   int  nAtom  = 0;    // Number of atoms on this processor.
   int  nGhost = 0;    // Number of ghosts on this processor.

   AtomIterator   atomIter;
   GhostIterator  ghostIter;

   double range = 0.02;
   //displaceAtoms(range);

   atomStorage.clearSnapshot();
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
   atomStorage.transformGenToCart(boundary);

   range = 0.02;
   for (int i=0; i < 4; ++i) {

      TEST_ASSERT(atomStorage.isCartesian());
      displaceAtoms(range);

      for (int j = 0; j < 2; ++j) {
         exchanger.update();
         TEST_ASSERT(nGhost == atomStorage.nGhost());
         TEST_ASSERT(nAtom == atomStorage.nAtom());
         displaceAtoms(range);
      }

      // Transform to generalized coordinates
      atomStorage.transformCartToGen(boundary);

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
      #ifdef SIMP_ANGLE
      if (hasAngles) {
         TEST_ASSERT(angleStorage.isValid(atomStorage, 
                     domain.communicator(), true));
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (hasDihedrals) {
         TEST_ASSERT(dihedralStorage.isValid(atomStorage, 
                     domain.communicator(), true));
      }
      #endif

      // Transform to Cartesian coordinates
      atomStorage.transformGenToCart(boundary);
      atomStorage.makeSnapshot();

   }

}

void ExchangerForceTest::testInitialForcesF() {
   printMethod(TEST_FUNC);
   reverseUpdateFlag = false;
   testInitialForces();
}

void ExchangerForceTest::testInitialForcesR() {
   printMethod(TEST_FUNC);
   reverseUpdateFlag = true;
   testInitialForces();
}

void ExchangerForceTest::testInitialForces()
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
   atomStorage.transformGenToCart(boundary);
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
   #ifdef SIMP_ANGLE
   if (hasAngles) {
      TEST_ASSERT(angleStorage.isValid(atomStorage, 
                  domain.communicator(), true));
   }
   #endif
   #ifdef SIMP_DIHEDRAL
   if (hasDihedrals) {
      TEST_ASSERT(dihedralStorage.isValid(atomStorage, 
                  domain.communicator(), true));
   }
   #endif

   TEST_ASSERT(pairPotential.reverseUpdateFlag() == reverseUpdateFlag);

   // Compute forces etc. with N^2 loop (MethodId = 2)
   pairPotential.setMethodId(2);
   computeForces();
   if (reverseUpdateFlag) {
      exchanger.reverseUpdate();
   }
   saveForces();
   int nPairNSq = 0;
   pairPotential.computeNPair(domain.communicator());
   if (domain.communicator().Get_rank() == 0) {
      nPairNSq = pairPotential.nPair();
   }
   double pairEnergyNSq = 0.0;
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
   int nPairList = 0;
   pairPotential.computeNPair(domain.communicator());
   if (domain.communicator().Get_rank() == 0) {
      nPairList = pairPotential.nPair();
   }
   double energyList = 0.0;
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
   Tensor bondStress;
   bondPotential.computeStress(domain.communicator());
   if (domain.communicator().Get_rank() == 0) {
      bondStress = bondPotential.stress();
   }
   #ifdef SIMP_ANGLE
   Tensor angleStress;
   if (hasDihedrals) {
      anglePotential.unsetStress();
      anglePotential.computeStress(domain.communicator());
      if (domain.communicator().Get_rank() == 0) {
         angleStress = anglePotential.stress();
      }
   }
   #endif
   #ifdef SIMP_DIHEDRAL
   Tensor dihedralStress;
   if (hasDihedrals) {
      dihedralPotential.unsetStress();
      dihedralPotential.computeStress(domain.communicator());
      if (domain.communicator().Get_rank() == 0) {
         dihedralStress = dihedralPotential.stress();
      }
   }
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
      for (i = 0; i < 3; ++i) {
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
   bondPotential.unsetStress();
   bondPotential.computeForcesAndStress(domain.communicator());
   #ifdef SIMP_ANGLE
   if (hasAngles) {
      anglePotential.unsetStress();
      anglePotential.computeForcesAndStress(domain.communicator());
   }
   #endif
   #ifdef SIMP_DIHEDRAL
   if (hasDihedrals) {
      dihedralPotential.unsetStress();
      dihedralPotential.computeForcesAndStress(domain.communicator());
   }
   #endif
   if (reverseUpdateFlag) {
      exchanger.reverseUpdate();
   }

   if (domain.communicator().Get_rank() == 0) {
      Tensor pairStress2 = pairPotential.stress();
      Tensor bondStress2;
      bondStress2 = bondPotential.stress();
      for (int i = 0; i < Dimension; ++i) {
         for (int j = 0; j < Dimension; ++j) {
            TEST_ASSERT(eq(pairStress(i, j), pairStress2(i, j)));
            TEST_ASSERT(eq(bondStress(i, j), bondStress2(i, j)));
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

void ExchangerForceTest::testForceCycleF() {
   printMethod(TEST_FUNC);
   reverseUpdateFlag = true;
   testForceCycle();
}

void ExchangerForceTest::testForceCycleR() {
   printMethod(TEST_FUNC);
   reverseUpdateFlag = false;
   testForceCycle();
}

void ExchangerForceTest::testForceCycle()
{
   initialize();

   // int  nAtom  = 0;    // Number of atoms on this processor.
   // int  nGhost = 0;    // Number of ghosts on this processor.
   bool needExchange;

   TEST_ASSERT(pairPotential.reverseUpdateFlag() == reverseUpdateFlag);

   // double range = 0.1;
   // displaceAtoms(range);

   atomStorage.clearSnapshot();
   TEST_ASSERT(!atomStorage.isCartesian());
   exchanger.exchange();
   exchangeNotify();

   // nAtom = atomStorage.nAtom();
   // nGhost = atomStorage.nGhost();

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
   atomStorage.transformGenToCart(boundary);
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
      needExchange = isExchangeNeeded(pairPotential.skin());
      //needExchange = atomStorage.needExchange(domain.communicator(), 
      //                                        pairPotential.skin());
      
      if (needExchange) {
         //if (domain.isMaster()) {
         //   std::cout << "Step " << i << ",  exchange " << j << std::endl;
         //}
         atomStorage.clearSnapshot();
         TEST_ASSERT(atomStorage.isCartesian());
         atomStorage.transformCartToGen(boundary);
         TEST_ASSERT(!atomStorage.isCartesian());

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

         // Build cell and pair lists
         pairPotential.buildCellList();
         TEST_ASSERT(!atomStorage.isCartesian());
         atomStorage.transformGenToCart(boundary);
         pairPotential.buildPairList();
         atomStorage.makeSnapshot();
         TEST_ASSERT(atomStorage.isCartesian());

         ++j;

      } else {

         exchanger.update();
      }

      // nAtom  = atomStorage.nAtom();
      // nGhost = atomStorage.nGhost();

      TEST_ASSERT(atomStorage.isValid());
      TEST_ASSERT(bondStorage.isValid(atomStorage, domain.communicator(),
                                      true)); 
      #ifdef SIMP_ANGLE
      if (hasAngles) {
         TEST_ASSERT(angleStorage.isValid(atomStorage, 
                     domain.communicator(), true));
      }
      #endif
      #ifdef SIMP_DIHEDRAL
      if (hasDihedrals) {
         TEST_ASSERT(dihedralStorage.isValid(atomStorage, 
                     domain.communicator(), true));
      }
      #endif

      // Calculate forces by an N^2 loop.
      computeForces();
      if (reverseUpdateFlag) {
         exchanger.reverseUpdate();
      }
      saveForces();
      nPairNSq = 0;
      pairEnergyNSq = 0.0;
      pairPotential.computeNPair(domain.communicator());
      pairPotential.computeEnergy(domain.communicator());
      if (domain.communicator().Get_rank() == 0) {
         nPairNSq = pairPotential.nPair();
         pairEnergyNSq = pairPotential.energy();
      }

      // Calculate forces etc. via pair list.   
      pairPotential.setMethodId(0); // PairList
      computeForces();
      if (reverseUpdateFlag) {
         exchanger.reverseUpdate();
      }
      energyList = 0.0;
      nPairList = 0;
      pairPotential.computeEnergy(domain.communicator());
      pairPotential.computeNPair(domain.communicator());
      if (domain.communicator().Get_rank() == 0) {
         energyList = pairPotential.energy();
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
         atomStorage.transformCartToGen(boundary);
         TEST_ASSERT(!atomStorage.isCartesian());
         pairPotential.buildCellList();
         atomStorage.transformGenToCart(boundary);
         pairPotential.buildPairList();
         TEST_ASSERT(atomStorage.isCartesian());
         pairPotential.setMethodId(0);    
         computeForces();
         saveForces();
         nPairF = 0;
         energyF = 0.0;
         pairPotential.computeEnergy(domain.communicator());
         pairPotential.computeNPair(domain.communicator());
         if (domain.communicator().Get_rank() == 0) {
            energyF = pairPotential.energy();
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
         atomStorage.transformGenToCart(boundary);
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

/*
* Determine whether an atom exchange and reneighboring is needed.
*/
bool ExchangerForceTest::isExchangeNeeded(double skin) 
{
   if (!atomStorage.isCartesian()) {
      UTIL_THROW("Error: Coordinates not Cartesian in isExchangeNeeded");
   } 

   // Calculate maximum square displacment on this node
   double maxSqDisp = atomStorage.maxSqDisplacement(); 
   int    needed = 0;
   if (sqrt(maxSqDisp) > 0.5*skin) {
      needed = 1; 
   }

   #if UTIL_MPI
   int neededAll;
   domain.communicator().Allreduce(&needed, &neededAll, 1, MPI::INT, MPI::MAX);
   return bool(neededAll);
   #else
   return bool(needed);
   #endif
}

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
