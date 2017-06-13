/*
* This program runs all unit tests in the mcMd directory.
*/

#include <test/CompositeTestRunner.h>

#include "chemistry/ChemistryTestComposite.h"
#include "neighbor/NeighborTestComposite.h"
#include "simulation/SimulationTestComposite.h"
#include "mcSimulation/McSimulationTest.h"
#include "mdSimulation/MdSimulationTest.h"

// Define a class McMdNsTestComposite
TEST_COMPOSITE_BEGIN(McMdNsTestComposite)
addChild(new ChemistryTestComposite, "chemistry/");
addChild(new NeighborTestComposite, "neighbor/");
addChild(new SimulationTestComposite, "simulation/");
addChild(new TEST_RUNNER(McSimulationTest), "mcSimulation/");
addChild(new TEST_RUNNER(MdSimulationTest), "mdSimulation/");
TEST_COMPOSITE_END

int main(int argc, char* argv[])
{
   if (argc > 2) {
      UTIL_THROW("Too many arguments");
   }

   McMdNsTestComposite runner;
   if (argc == 2) {
      runner.addFilePrefix(argv[1]);
    }
   runner.run();

}
