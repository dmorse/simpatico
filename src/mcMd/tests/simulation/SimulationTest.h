#ifndef MCMD_SIMULATION_TEST_H
#define MCMD_SIMULATION_TEST_H

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
      openFile("in/SimulationAngleDihedral"); 
      #else
      openFile("in/SimulationAngle"); 
      #endif
      #else
      #ifdef SIMP_DIHEDRAL
      openFile("in/SimulationDihedral"); 
      #else
      openFile("in/Simulation"); 
      #endif
      #endif
      simulation_.readParameters(file());
   }

   void testReadParam();

private:

   Simulation simulation_;

};



void SimulationTest::testReadParam()
{
   printMethod(TEST_FUNC);

   if (verbose() > 1) {
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
