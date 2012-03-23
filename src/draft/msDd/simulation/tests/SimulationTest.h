#ifndef MSDD_SIMULATION_TEST_H
#define MSDD_SIMULATION_TEST_H

#include <msDd/simulation/Simulation.h>

#ifdef UTIL_MPI
#define TEST_MPI
#endif

#include <test/ParamFileTest.h>
#include <test/UnitTestRunner.h>

using namespace Util;

class SimulationTest : public UnitTest
{

   MsDd::Simulation simulation;

public:

   SimulationTest()
    : simulation(MPI::COMM_WORLD)
   {}

   virtual void setUp()
   {}

   void testReadParam();

};

inline void SimulationTest::testReadParam()
{  
   printMethod(TEST_FUNC); 

   std::ifstream paramFile;
   if (mpiRank() == 0) {
      openInputFile("in/param", paramFile); 
      simulation.readParam(paramFile); 
      //if (verbose() > 0) {
      simulation.writeParam(std::cout);
      //}
   }
   // TEST_ASSERT(simulation.buffer().isAllocated())
   //TEST_ASSERT(simulation.domain().isInitialized());
   //TEST_ASSERT(simulation.domain().hasBoundary());
}

TEST_BEGIN(SimulationTest)
TEST_ADD(SimulationTest, testReadParam)
TEST_END(SimulationTest)

#endif
