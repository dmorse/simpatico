#ifndef NEIGHBOR_TEST_COMPOSITE_H
#define NEIGHBOR_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "CellTest.h"
#include "CellListTest.h"
#include "PairListTest.h"

TEST_COMPOSITE_BEGIN(NeighborTestComposite)
TEST_COMPOSITE_ADD_UNIT(CellTest);
TEST_COMPOSITE_ADD_UNIT(CellListTest);
TEST_COMPOSITE_ADD_UNIT(PairListTest);
TEST_COMPOSITE_END

#endif
