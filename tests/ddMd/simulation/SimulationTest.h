#ifndef DDMD_SYSTEM_TEST_H
#define DDMD_SYSTEM_TEST_H

#include <ddMd/simulation/Simulation.h>
#include <ddMd/storage/AtomIterator.h>
#include <ddMd/storage/GhostIterator.h>
#include <ddMd/potentials/pair/PairPotential.h>
#include <ddMd/integrators/Integrator.h>
#include <util/random/Random.h>
#include <util/format/Dbl.h>
#include <util/mpi/MpiLogger.h>

#ifdef UTIL_MPI
#define TEST_MPI
#endif

#include <test/ParamFileTest.h>
#include <test/UnitTestRunner.h>

using namespace Util;
using namespace DdMd;

class SimulationTest : public ParamFileTest<Simulation>
{

public:

   virtual void setUp()
   {}

   void displaceAtoms(AtomStorage& atomStorage, const Boundary& boundary, 
                      Random& random, double range);

   void testReadParam();

   void testReadConfig();

   void testExchangeAtoms();

   void testExchange();

   void testUpdate();

   void testCalculateForces();

   void testIntegrate1();

};


inline void SimulationTest::displaceAtoms(AtomStorage& atomStorage, 
                                          const Boundary& boundary, 
                                          Random& random, double range)
{
   Vector ranges;
   double min, max;
   if (UTIL_ORTHOGONAL) {
      for (int i = 0; i < Dimension; ++i) {
         ranges[i] = range;
      }
   } else {
      for (int i = 0; i < Dimension; ++i) {
         ranges[i] = range/boundary.length(i);
      }
   }
   AtomIterator atomIter;
   for(int i = 0; i < Dimension; ++i) {
      max = ranges[i];
      min = -max;
      atomStorage.begin(atomIter);
      for ( ; atomIter.notEnd(); ++atomIter) {
         atomIter->position()[i] += random.uniform(min, max);
      }
   }
}

inline void SimulationTest::testReadParam()
{  
   printMethod(TEST_FUNC); 

   openFile("in/param1"); 
   object().readParam(file()); 
   if (verbose() > 0) {
      object().writeParam(std::cout);
   }
   // TEST_ASSERT(object().buffer().isAllocated())
   TEST_ASSERT(object().domain().isInitialized());
   TEST_ASSERT(object().domain().hasBoundary());
}

inline void SimulationTest::testReadConfig()
{
   printMethod(TEST_FUNC); 

   openFile("in/param1"); 
   object().readParam(file()); 

   Domain& domain = object().domain();
   Boundary& boundary = object().boundary();
   AtomStorage& atomStorage = object().atomStorage();

   std::string filename("config1");
   object().readConfig(filename);

   int nAtomAll;
   int myRank = domain.gridRank();
   atomStorage.computeNAtomTotal(domain.communicator());
   if (myRank == 0) {
      nAtomAll = atomStorage.nAtomTotal();
      //std::cout << "Total atom count = " << nAtomAll << std::endl;
      TEST_ASSERT(nAtomAll == 100);
   }
   std::cout.flush();

   // Check that all atoms are within the processor domain.
   int j = 0;
   AtomIterator atomIter;
   atomStorage.begin(atomIter);
   for ( ; atomIter.notEnd(); ++atomIter) {
      j++;
      TEST_ASSERT( domain.isInDomain( atomIter->position() ) );
   }
   TEST_ASSERT(j == atomStorage.nAtom());

}

inline void SimulationTest::testExchangeAtoms()
{
   printMethod(TEST_FUNC); 

   openFile("in/param1"); 
   object().readParam(file()); 

   Domain&   domain = object().domain();
   Boundary& boundary = object().boundary();
   AtomStorage& atomStorage = object().atomStorage();
   Random&  random = object().random();

   std::string filename("config1");
   object().readConfig(filename);

   displaceAtoms(atomStorage, boundary, random, 0.8);
 
   #if 0
   // Range of random increments for the atom positions.
   double range1 = double(-0.8);
   double range2 = double(0.8);

   // Add a random increment to atom positions
   AtomIterator atomIter;
   for (int i = 0; i < Dimension; ++i) {
     atomStorage.begin(atomIter);
     for ( ; atomIter.notEnd(); ++atomIter) {
        atomIter->position()[i] += random.uniform(range1, range2);
     }
   }
   #endif

   // Exchange atoms among processors
   object().exchanger().exchange();

   // Check that all atoms are accounted for after exchange.
   int nAtomAll;
   int myRank   = domain.gridRank();
   atomStorage.computeNAtomTotal(domain.communicator());
   if (myRank == 0) {
      nAtomAll = atomStorage.nAtomTotal();
      //std::cout << "Total atom count = " << nAtomAll << std::endl;
      TEST_ASSERT(nAtomAll == 100);
   }

   // Check that all atoms are within the processor domain.
   int j = 0;
   AtomIterator atomIter;
   atomStorage.begin(atomIter);
   for ( ; atomIter.notEnd(); ++atomIter) {
      j++;
      TEST_ASSERT( domain.isInDomain( atomIter->position() ) );
   }
   TEST_ASSERT(j == atomStorage.nAtom());

}

inline void SimulationTest::testExchange()
{
   printMethod(TEST_FUNC); 

   openFile("in/param1"); 
   object().readParam(file()); 

   Domain&  domain  = object().domain();
   Boundary& boundary = object().boundary();
   AtomStorage& atomStorage = object().atomStorage();
   Random&  random = object().random();

   std::string filename("config1");
   object().readConfig(filename);

   displaceAtoms(atomStorage, boundary, random, 0.8);

   // Exchange atoms among processors
   object().exchanger().exchange();

   // Check that all atoms are accounted for after exchange.
   int nAtomAll;
   int myRank = domain.gridRank();
   atomStorage.computeNAtomTotal(domain.communicator());
   if (myRank == 0) {
      nAtomAll = atomStorage.nAtomTotal();
      //std::cout << "Total atom count = " << nAtomAll << std::endl;
      TEST_ASSERT(nAtomAll == 100);
   }

   // Check that all atoms are within the processor domain.
   AtomIterator atomIter;
   int j = 0;
   atomStorage.begin(atomIter);
   for ( ; atomIter.notEnd(); ++atomIter) {
      j++;
      TEST_ASSERT(domain.isInDomain( atomIter->position()));
   }
   TEST_ASSERT(j == atomStorage.nAtom());
   int nAtom = atomStorage.nAtom();

   // Check that all ghosts are outside the processor domain.
   GhostIterator ghostIter;
   atomStorage.begin(ghostIter);
   for ( ; ghostIter.notEnd(); ++ghostIter) {
      TEST_ASSERT(!domain.isInDomain(ghostIter->position()));
   }

   #if 0
   int nGhostAll = object().nGhostTotal();
   if (myRank == 0) {
      std::cout << "Total ghost count = " << nGhostAll << std::endl;
   }
   #endif

}

inline void SimulationTest::testUpdate()
{
   printMethod(TEST_FUNC); 

   openFile("in/param1"); 
   object().readParam(file()); 

   Domain&  domain  = object().domain();
   Boundary& boundary = object().boundary();
   AtomStorage& atomStorage = object().atomStorage();
   Random&  random  = object().random();

   std::string filename("config1");
   object().readConfig(filename);

   displaceAtoms(atomStorage, boundary, random, 0.8);

   // Exchange atoms among processors
   object().exchanger().exchange();

   // Check that all atoms are accounted for after exchange.
   int nAtomAll;
   int myRank = domain.gridRank();
   atomStorage.computeNAtomTotal(domain.communicator());
   if (myRank == 0) {
      nAtomAll = atomStorage.nAtomTotal();
      //std::cout << "Total atom count = " << nAtomAll << std::endl;
      TEST_ASSERT(nAtomAll == 100);
   }

   // Check that all atoms are within the processor domain.
   int j = 0;
   AtomIterator atomIter;
   atomStorage.begin(atomIter);
   for ( ; atomIter.notEnd(); ++atomIter) {
      j++;
      TEST_ASSERT( domain.isInDomain( atomIter->position() ) );
   }
   TEST_ASSERT(j == atomStorage.nAtom());
   int nAtom = atomStorage.nAtom();

   // Check that all ghosts are outside the processor domain.
   GhostIterator ghostIter;
   atomStorage.begin(ghostIter);
   for ( ; ghostIter.notEnd(); ++ghostIter) {
      TEST_ASSERT(!domain.isInDomain(ghostIter->position()));
   }

   if (!UTIL_ORTHOGONAL) {
       atomStorage.transformGenToCart(boundary);
   }

   displaceAtoms(atomStorage, boundary, random, 0.1);
   object().exchanger().update();

}

inline void SimulationTest::testCalculateForces()
{
   printMethod(TEST_FUNC); 

   openFile("in/param1"); 
   object().readParam(file()); 

   Domain&  domain  = object().domain();
   Boundary& boundary = object().boundary();
   AtomStorage& atomStorage = object().atomStorage();
   int myRank = domain.gridRank();

   std::string filename("config1");
   object().readConfig(filename);

   // Compute forces.
   object().pairPotential().buildCellList();
   if (!UTIL_ORTHOGONAL) {
      atomStorage.transformGenToCart(object().boundary());
   }
   object().pairPotential().buildPairList();
   object().computeForces();

   Vector f(0.0); // total force on processor
   Vector t(0.0); // total force on all processors
   AtomIterator iter;
   atomStorage.begin(iter);
   for ( ; iter.notEnd(); ++iter) {
      f += iter->force();
   }

   #ifdef UTIL_MPI
   domain.communicator().Reduce(&f[0], &t[0], 1, 
                                 MPI::DOUBLE, MPI::SUM, 0);
   domain.communicator().Reduce(&f[1], &t[1], 1, 
                                 MPI::DOUBLE, MPI::SUM, 0);
   domain.communicator().Reduce(&f[2], &t[2], 1, 
                                 MPI::DOUBLE, MPI::SUM, 0);

   if (myRank == 0) {
     std::cout << std::endl;
     std::cout << t << std::endl;
   }
   #endif

}

inline void SimulationTest::testIntegrate1()
{
   printMethod(TEST_FUNC); 

   openFile("in/param2"); 
   object().readParam(file()); 

   Domain& domain = object().domain();
   Boundary& boundary = object().boundary();
   AtomStorage& atomStorage = object().atomStorage();
   int myRank = domain.gridRank();

   // Read configuration file
   std::string filename("config2");
   object().readConfig(filename);

   // Set random velocities
   double temperature = 1.0;
   object().setBoltzmannVelocities(temperature);

   if (myRank == 0) {
      std::cout << std::endl;
   }

   // Setup the integrator
   object().integrator().setup();
   TEST_ASSERT(object().isValid());

   double kinetic;
   double potential;
   for (int i = 0; i < 10; ++i ) {

      // Calculate energies
      object().computeKineticEnergy();
      object().computePotentialEnergies();
      if (myRank == 0) {
         kinetic = object().kineticEnergy();
         potential = object().potentialEnergy();
         std::cout << Dbl(kinetic) << Dbl(potential) 
                   << Dbl(kinetic + potential) << std::endl;
      }

      object().integrator().run(500);
      TEST_ASSERT(object().isValid());
   }
   object().integrator().outputStatistics(std::cout);

}

TEST_BEGIN(SimulationTest)
TEST_ADD(SimulationTest, testReadParam)
TEST_ADD(SimulationTest, testReadConfig)
TEST_ADD(SimulationTest, testExchangeAtoms)
TEST_ADD(SimulationTest, testExchange)
TEST_ADD(SimulationTest, testUpdate)
TEST_ADD(SimulationTest, testCalculateForces)
TEST_ADD(SimulationTest, testIntegrate1)
TEST_END(SimulationTest)

#endif
