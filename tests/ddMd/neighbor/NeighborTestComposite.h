#ifndef DDMD_NEIGHBOR_TEST_COMPOSITE_H
#define DDMD_NEIGHBOR_TEST_COMPOSITE_H

#include "CellTest.h"
#include "CellListTest.h"
#include "PairListTest.h"

#include <test/CompositeTestRunner.h>

TEST_COMPOSITE_BEGIN(NeighborTestComposite)
TEST_COMPOSITE_ADD_UNIT(CellTest);
TEST_COMPOSITE_ADD_UNIT(CellListTest);
TEST_COMPOSITE_ADD_UNIT(PairListTest);
TEST_COMPOSITE_END

#endif
