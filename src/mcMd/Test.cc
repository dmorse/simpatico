#ifndef MCMD_TEST_CC
#define MCMD_TEST_CC

/*
* This program runs all unit tests in the mcMd directory.
*/

#include <test/CompositeTestRunner.h>

#include <mcMd/boundary/tests/BoundaryTestComposite.h>
#include <mcMd/chemistry/tests/ChemistryTestComposite.h>
#include <mcMd/neighbor/tests/NeighborTestComposite.h>
#include <mcMd/simulation/tests/SimulationTestComposite.h>
#include <mcMd/potentials/tests/PotentialTestComposite.h>
#include <mcMd/mcSimulation/tests/McSimulationTest.h>
#include <mcMd/mdSimulation/tests/MdSimulationTest.h>

// Define a class McMdNsTestComposite
TEST_COMPOSITE_BEGIN(McMdNsTestComposite)
addChild(new BoundaryTestComposite, "boundary/tests/");
addChild(new ChemistryTestComposite, "chemistry/tests/");
addChild(new NeighborTestComposite, "neighbor/tests/");
addChild(new SimulationTestComposite, "simulation/tests/");
addChild(new PotentialTestComposite, "potentials/tests/");
addChild(new TEST_RUNNER(McSimulationTest), "mcSimulation/tests/");
addChild(new TEST_RUNNER(MdSimulationTest), "mdSimulation/tests/");
TEST_COMPOSITE_END

int main() {

   McMdNsTestComposite runner;
   runner.run();

}
#endif
