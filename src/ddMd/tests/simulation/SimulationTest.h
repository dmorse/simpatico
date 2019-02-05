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
#include <util/misc/FileMaster.h>

#ifdef UTIL_MPI
#ifndef TEST_MPI
#define TEST_MPI
#endif
#endif

#include <test/ParamFileTest.h>
#include <test/UnitTestRunner.h>
#include <test/CommandLine.h>

using namespace Util;
using namespace DdMd;

class SimulationTest : public ParamFileTest
{
private:

   DdMd::Simulation simulation_;

   void displaceAtoms(AtomStorage& atomStorage, const Boundary& boundary, 
                      Random& random, double range);

public:

   virtual void setUp()
   {  
      Label::clear();
      simulation_.fileMaster().setRootPrefix(filePrefix()); 
   }

   virtual void tearDown()
   {  
      Label::clear();
   }

   void testReadParam();

   void testReadConfig();

   void testExchangeAtoms();

   void testExchange();

   void testUpdate();

   void testCalculateForces();

   void testIntegrate1();

};


inline void 
SimulationTest::displaceAtoms(AtomStorage& atomStorage, 
                              const Boundary& boundary, 
                              Random& random, double range)
{
   Vector ranges;
   double min, max;
   for (int i = 0; i < Dimension; ++i) {
      ranges[i] = range/boundary.length(i);
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
   simulation_.readParam(file()); 
   if (verbose() > 0) {
      simulation_.writeParam(std::cout);
   }
   // TEST_ASSERT(simulation_.buffer().isAllocated())
   TEST_ASSERT(simulation_.domain().isInitialized());
   TEST_ASSERT(simulation_.domain().hasBoundary());
}

inline void SimulationTest::testReadConfig()
{
   printMethod(TEST_FUNC); 

   openFile("in/param1"); 
   simulation_.readParam(file()); 

   Domain& domain = simulation_.domain();
   //Boundary& boundary = simulation_.boundary();
   AtomStorage& atomStorage = simulation_.atomStorage();

   std::string filename("config1");
   simulation_.readConfig(filename);

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
   simulation_.readParam(file()); 
   file().close(); 

   Domain&   domain = simulation_.domain();
   Boundary& boundary = simulation_.boundary();
   AtomStorage& atomStorage = simulation_.atomStorage();
   Random&  random = simulation_.random();

   std::string filename("config1");
   simulation_.readConfig(filename);

   displaceAtoms(atomStorage, boundary, random, 0.8);
 
   // Exchange atoms among processors
   simulation_.exchanger().exchange();

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
      TEST_ASSERT( domain.isInDomain(atomIter->position()) );
   }
   TEST_ASSERT(j == atomStorage.nAtom());

}

inline void SimulationTest::testExchange()
{
   printMethod(TEST_FUNC); 

   openFile("in/param1"); 
   simulation_.readParam(file());
   file().close(); 

   Boundary& boundary = simulation_.boundary();
   Random& random = simulation_.random();
   AtomStorage& atomStorage = simulation_.atomStorage();

   std::string filename("config1");
   simulation_.readConfig(filename);

   displaceAtoms(atomStorage, boundary, random, 0.8);

   // Exchange atoms among processors
   simulation_.exchanger().exchange();
}

inline void SimulationTest::testUpdate()
{
   printMethod(TEST_FUNC); 

   openFile("in/param1"); 
   simulation_.readParam(file()); 

   Domain&  domain  = simulation_.domain();
   Boundary& boundary = simulation_.boundary();
   AtomStorage& atomStorage = simulation_.atomStorage();
   Random&  random  = simulation_.random();

   std::string filename("config1");
   simulation_.readConfig(filename);

   displaceAtoms(atomStorage, boundary, random, 0.8);

   // Exchange atoms among processors
   simulation_.exchanger().exchange();

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

   // Check that all ghosts are outside the processor domain.
   GhostIterator ghostIter;
   atomStorage.begin(ghostIter);
   for ( ; ghostIter.notEnd(); ++ghostIter) {
      TEST_ASSERT(!domain.isInDomain(ghostIter->position()));
   }
   atomStorage.transformGenToCart(boundary);

   #if 1
   displaceAtoms(atomStorage, boundary, random, 0.1);
   simulation_.exchanger().update();
   #endif

}

inline void SimulationTest::testCalculateForces()
{
   printMethod(TEST_FUNC); 

   openFile("in/param1"); 
   simulation_.readParam(file()); 

   Domain&  domain  = simulation_.domain();
   AtomStorage& atomStorage = simulation_.atomStorage();
   int myRank = domain.gridRank();

   std::string filename("config1");
   simulation_.readConfig(filename);

   // Compute forces.
   simulation_.pairPotential().buildCellList();
   atomStorage.transformGenToCart(simulation_.boundary());
   simulation_.pairPotential().buildPairList();
   simulation_.computeForces();

   Vector f(0.0); // total force on processor
   Vector t(0.0); // total force on all processors
   AtomIterator iter;
   atomStorage.begin(iter);
   for ( ; iter.notEnd(); ++iter) {
      f += iter->force();
   }

   #ifdef UTIL_MPI
   domain.communicator().Reduce(&f[0], &t[0], 1, 
                                 MPI_DOUBLE, MPI_SUM, 0);
   domain.communicator().Reduce(&f[1], &t[1], 1, 
                                 MPI_DOUBLE, MPI_SUM, 0);
   domain.communicator().Reduce(&f[2], &t[2], 1, 
                                 MPI_DOUBLE, MPI_SUM, 0);

   if (myRank == 0) {
     std::cout << std::endl;
     std::cout << t << std::endl;
   }
   #endif

}

inline void SimulationTest::testIntegrate1()
{
   printMethod(TEST_FUNC); 

   Domain& domain = simulation_.domain();

   CommandLine opts;
   opts.append("-e");
   simulation_.setOptions(opts.argc(), opts.argv());

   openFile("in/param2"); 
   simulation_.readParam(file()); 

   int myRank = domain.gridRank();
   if (myRank == 0) {
      std::cout << std::endl;
   }

   // Read configuration file
   std::string filename("config2");
   simulation_.readConfig(filename);

   // Set random velocities
   double temperature = 1.0;
   simulation_.setBoltzmannVelocities(temperature);

   TEST_ASSERT(simulation_.isValid());

   double kinetic;
   double potential;
   for (int i = 0; i < 4; ++i ) {

      // Calculate energies
      simulation_.computeKineticEnergy();
      simulation_.computePotentialEnergies();
      if (myRank == 0) {
         kinetic = simulation_.kineticEnergy();
         potential = simulation_.potentialEnergy();
         std::cout << Dbl(kinetic) << Dbl(potential) 
                   << Dbl(kinetic + potential) << std::endl;
      }

      simulation_.integrator().run(20);
      TEST_ASSERT(simulation_.isValid());
   }
   if (myRank == 0) {
      simulation_.integrator().outputStatistics(std::cout);
   }

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
