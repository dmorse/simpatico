#ifndef MCMD_SIMULATION_TEST_COMPOSITE_H
#define MCMD_SIMULATION_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "SystemTest.h"
#include "SimulationTest.h"

TEST_COMPOSITE_BEGIN(SimulationTestComposite)
TEST_COMPOSITE_ADD_UNIT(SystemTest)
TEST_COMPOSITE_ADD_UNIT(SimulationTest)
TEST_COMPOSITE_END

#endif
