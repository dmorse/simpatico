#ifndef DDMD_SP_NEIGHBOR_TEST_COMPOSITE_H
#define DDMD_SP_NEIGHBOR_TEST_COMPOSITE_H

#include "SpCellTest.h"
#include "SpCellListTest.h"

#include <test/CompositeTestRunner.h>

TEST_COMPOSITE_BEGIN(SpNeighborTestComposite)
TEST_COMPOSITE_ADD_UNIT(SpCellTest);
TEST_COMPOSITE_ADD_UNIT(SpCellListTest);
TEST_COMPOSITE_END

#endif
