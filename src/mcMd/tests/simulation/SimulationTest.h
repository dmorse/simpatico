#ifndef MCMD_SIMULATION_TEST_H
#define MCMD_SIMULATION_TEST_H

#ifdef UTIL_MPI
#define TEST_MPI
#endif

#include <test/ParamFileTest.h>
#include <test/UnitTestRunner.h>

#include <mcMd/simulation/Simulation.h>
#include <simp/species/Species.h>
#include <mcMd/chemistry/Molecule.h>

using namespace Util;
using namespace McMd;

class SimulationTest : public ParamFileTest
{

public:

   SimulationTest()
   {}

   virtual void setUp()
   {  
      #ifdef SIMP_ANGLE
      #ifdef SIMP_DIHEDRAL
      openInputFile("in/SimulationAngleDihedral", file()); 
      #else
      openInputFile("in/SimulationAngle", file()); 
      #endif
      #else
      #ifdef SIMP_DIHEDRAL
      openInputFile("in/SimulationDihedral", file()); 
      #else
      openInputFile("in/Simulation", file()); 
      #endif
      #endif
   }

   void testReadParam();

private:

   Simulation simulation_;

};



void SimulationTest::testReadParam()
{
   printMethod(TEST_FUNC);

   simulation_.readParameters(file());
   if (verbose() > 1 && isIoProcessor()) {
      std::cout << std::endl;
      simulation_.writeParam(std::cout);
   }

   Species *species;
   for (int i=0; i < simulation_.nSpecies() ; ++i ) {

      species = &simulation_.species(i);
      TEST_ASSERT(species->isValid());

   }

   TEST_ASSERT(simulation_.isValid());
}

TEST_BEGIN(SimulationTest)
TEST_ADD(SimulationTest, testReadParam)
TEST_END(SimulationTest)

#endif
