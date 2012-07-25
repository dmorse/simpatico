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

   void displaceAtoms(AtomStorage& storage, const Boundary& boundary, Random& random, double range);

   void testReadParam();

   void testReadConfig();

   void testExchangeAtoms();

   void testExchange();

   void testUpdate();

   void testCalculateForces();

   void testIntegrate1();

};


inline void SimulationTest::displaceAtoms(AtomStorage& storage, const Boundary& boundary, 
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
      storage.begin(atomIter);
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
   AtomStorage& storage = object().atomStorage();

   std::string filename("config1");
   object().readConfig(filename);

   int nAtomAll;
   int myRank = domain.gridRank();
   storage.computeNAtomTotal(domain.communicator());
   if (myRank == 0) {
      nAtomAll = storage.nAtomTotal();
      //std::cout << "Total atom count = " << nAtomAll << std::endl;
      TEST_ASSERT(nAtomAll == 100);
   }
   std::cout.flush();

   #if 0
   MpiLogger logger;
   logger.begin();

   std::cout << std::endl;
   std::cout << "Processor " << myRank << std::endl;
   std::cout << "Position  ";
   for (int i = 0; i < Dimension; ++i) {
       std::cout << domain.gridCoordinate(i) << "  ";
   }
   std::cout << std::endl;

   for (int i = 0; i < Dimension; ++i) {
      std::cout << "bound(" << i << ", 0) = " 
                << Dbl(domain.domainBound(i, 0));
      std::cout << "  bound(" << i << ", 1) = " 
                << Dbl(domain.domainBound(i, 1));
      std::cout << std::endl;
   }
   std::cout << std::endl;
   std::cout << "Lengths " << object().boundary().lengths()
             << std::endl;
   std::cout << "nAtom  " << storage.nAtom() << std::endl;
   std::cout << "nGhost " << storage.nGhost()<< std::endl;
   std::cout << std::endl;
   logger.end();
   #endif

   // Check that all atoms are within the processor domain.
   int j = 0;
   AtomIterator atomIter;
   storage.begin(atomIter);
   for ( ; atomIter.notEnd(); ++atomIter) {
      j++;
      TEST_ASSERT( domain.isInDomain( atomIter->position() ) );
   }
   TEST_ASSERT(j == storage.nAtom());

}

inline void SimulationTest::testExchangeAtoms()
{
   printMethod(TEST_FUNC); 

   openFile("in/param1"); 
   object().readParam(file()); 

   Domain&   domain = object().domain();
   Boundary& boundary = object().boundary();
   AtomStorage& storage = object().atomStorage();
   Random&  random = object().random();

   std::string filename("config1");
   object().readConfig(filename);

   displaceAtoms(storage, boundary, random, 0.8);
 
   #if 0
   // Range of random increments for the atom positions.
   double range1 = double(-0.8);
   double range2 = double(0.8);

   // Add a random increment to atom positions
   AtomIterator atomIter;
   for (int i = 0; i < Dimension; ++i) {
     storage.begin(atomIter);
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
   storage.computeNAtomTotal(domain.communicator());
   if (myRank == 0) {
      nAtomAll = storage.nAtomTotal();
      //std::cout << "Total atom count = " << nAtomAll << std::endl;
      TEST_ASSERT(nAtomAll == 100);
   }

   // Check that all atoms are within the processor domain.
   int j = 0;
   AtomIterator atomIter;
   storage.begin(atomIter);
   for ( ; atomIter.notEnd(); ++atomIter) {
      j++;
      TEST_ASSERT( domain.isInDomain( atomIter->position() ) );
   }
   TEST_ASSERT(j == storage.nAtom());

}

inline void SimulationTest::testExchange()
{
   printMethod(TEST_FUNC); 

   openFile("in/param1"); 
   object().readParam(file()); 

   Domain&  domain  = object().domain();
   Boundary& boundary = object().boundary();
   AtomStorage& storage = object().atomStorage();
   Random&  random = object().random();

   std::string filename("config1");
   object().readConfig(filename);

   displaceAtoms(storage, boundary, random, 0.8);

   #if 0
   // Range of random increments for the atom positions.
   double range1 = double(-0.8);
   double range2 = double(0.8);

   // Add a random increment to atom positions
   AtomIterator atomIter;
   for (int i = 0; i < Dimension; ++i) {
     storage.begin(atomIter);
     for ( ; atomIter.notEnd(); ++atomIter) {
        atomIter->position()[i] += random.uniform(range1, range2);
     }
   }
   #endif

   // Exchange atoms among processors
   object().exchanger().exchange();

   // Check that all atoms are accounted for after exchange.
   int nAtomAll;
   int myRank = domain.gridRank();
   storage.computeNAtomTotal(domain.communicator());
   if (myRank == 0) {
      nAtomAll = storage.nAtomTotal();
      //std::cout << "Total atom count = " << nAtomAll << std::endl;
      TEST_ASSERT(nAtomAll == 100);
   }

   // Check that all atoms are within the processor domain.
   AtomIterator atomIter;
   int j = 0;
   storage.begin(atomIter);
   for ( ; atomIter.notEnd(); ++atomIter) {
      j++;
      TEST_ASSERT(domain.isInDomain( atomIter->position()));
   }
   TEST_ASSERT(j == storage.nAtom());
   int nAtom = storage.nAtom();

   // Check that all ghosts are outside the processor domain.
   GhostIterator ghostIter;
   storage.begin(ghostIter);
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
   AtomStorage& storage = object().atomStorage();
   Random&  random  = object().random();

   std::string filename("config1");
   object().readConfig(filename);

   displaceAtoms(storage, boundary, random, 0.8);

   #if 0
   // Add a random increment to atom positions
   double range1 = double(-0.8);
   double range2 = double(0.8);
   AtomIterator atomIter;
   for (int i = 0; i < Dimension; ++i) {
     storage.begin(atomIter);
     for ( ; atomIter.notEnd(); ++atomIter) {
        atomIter->position()[i] += random.uniform(range1, range2);
     }
   }
   #endif

   // Exchange atoms among processors
   object().exchanger().exchange();

   // Check that all atoms are accounted for after exchange.
   int nAtomAll;
   int myRank = domain.gridRank();
   storage.computeNAtomTotal(domain.communicator());
   if (myRank == 0) {
      nAtomAll = storage.nAtomTotal();
      //std::cout << "Total atom count = " << nAtomAll << std::endl;
      TEST_ASSERT(nAtomAll == 100);
   }


   // Check that all atoms are within the processor domain.
   int j = 0;
   AtomIterator atomIter;
   storage.begin(atomIter);
   for ( ; atomIter.notEnd(); ++atomIter) {
      j++;
      TEST_ASSERT( domain.isInDomain( atomIter->position() ) );
   }
   TEST_ASSERT(j == storage.nAtom());
   int nAtom = storage.nAtom();

   // Check that all ghosts are outside the processor domain.
   GhostIterator ghostIter;
   storage.begin(ghostIter);
   for ( ; ghostIter.notEnd(); ++ghostIter) {
      TEST_ASSERT(!domain.isInDomain(ghostIter->position()));
   }

   #if 0
   int nGhostAll = object().nGhostTotal();
   if (myRank == 0) {
      nGhostAll = object().nGhostTotal();
      std::cout << "Total ghost count = " << nGhostAll << std::endl;
   }
   #endif

   #if 0
   //Print number of atoms on each processor after the ghost exchange.
   MpiLogger logger;
   logger.begin();
   std::cout << "Processor " << myRank 
             << " : Post-ghost exchange Atoms count = "
             << storage.nAtom() << std::endl;
   logger.end();

   // Print number of ghosts on each processor after the exchange.
   logger.begin();
   std::cout << "Processor " << myRank 
             << " : Post-ghost exchange Ghost count = "
             << storage.nGhost() << std::endl;
   logger.end();
   #endif

   displaceAtoms(storage, boundary, random, 0.1);

   #if 0
   // Add a random increment to atom positions
   range1 = double(-0.1);
   range2 = double(+0.1);
   storage.begin(atomIter);
   for ( ; atomIter.notEnd(); ++atomIter) {
      for (int i = 0; i < Dimension; ++i) {
         atomIter->position()[i] += random.uniform(range1, range2);
      }
   }
   #endif

   object().exchanger().update();

}

inline void SimulationTest::testCalculateForces()
{
   printMethod(TEST_FUNC); 

   openFile("in/param1"); 
   object().readParam(file()); 

   Domain&  domain  = object().domain();
   Boundary& boundary = object().boundary();
   AtomStorage& storage = object().atomStorage();
   int myRank = domain.gridRank();

   std::string filename("config1");
   object().readConfig(filename);
   //object().exchanger().exchange();

   // Compute forces.
   object().pairPotential().buildCellList();
   if (!UTIL_ORTHOGONAL) {
      storage.transformGenToCart(object().boundary());
   }
   object().pairPotential().buildPairList();
   object().computeForces();

   Vector f(0.0); // total force on processor
   Vector t(0.0); // total force on all processors
   AtomIterator iter;
   storage.begin(iter);
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
   AtomStorage& storage = object().atomStorage();
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
