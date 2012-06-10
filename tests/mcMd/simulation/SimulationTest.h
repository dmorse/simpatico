#ifndef MCMD_SIMULATION_TEST_H
#define MCMD_SIMULATION_TEST_H

#include <test/ParamFileTest.h>
#include <test/UnitTestRunner.h>

#include <mcMd/simulation/Simulation.h>
#include <mcMd/species/Species.h>
#include <mcMd/chemistry/Molecule.h>

using namespace Util;
using namespace McMd;

class SimulationTest : public ParamFileTest<Simulation>
{

public:

   SimulationTest()
   {}

   virtual void setUp()
   {  

      #ifdef INTER_ANGLE
      #ifdef INTER_DIHEDRAL
      openFile("in/SimulationAngleDihedral"); 
      #else
      openFile("in/SimulationAngle"); 
      #endif
      #else
      #ifdef INTER_DIHEDRAL
      openFile("in/SimulationDihedral"); 
      #else
      openFile("in/Simulation"); 
      #endif
      #endif
      object().readParam(file());
   }

   void testReadParam();

};



void SimulationTest::testReadParam()
{
   printMethod(TEST_FUNC);

   if (verbose() > 1) {
      std::cout << std::endl;
      object().writeParam(std::cout);
   }

   Species *species;
   for (int i=0; i < object().nSpecies() ; ++i ) {

      species = &object().species(i);
      TEST_ASSERT(species->isValid());

   }

   TEST_ASSERT(object().isValid());
}

TEST_BEGIN(SimulationTest)
TEST_ADD(SimulationTest, testReadParam)
TEST_END(SimulationTest)

#endif
