#ifndef DDMD_MODIFIER_TEST_COMPOSITE_H
#define DDMD_MODIFIER_TEST_COMPOSITE_H

#include "ModifierTest.h"

#include <test/CompositeTestRunner.h>

TEST_COMPOSITE_BEGIN(NeighborTestComposite)
TEST_COMPOSITE_ADD_UNIT(ModifierTest);
TEST_COMPOSITE_END

#endif
