#ifndef DDMD_TESTS_MODIFIER_TEST_COMPOSITE_H
#define DDMD_TESTS_MODIFIER_TEST_COMPOSITE_H

#include "ModifierTest.h"
#include "ModifierManagerTest.h"

#include <test/CompositeTestRunner.h>

TEST_COMPOSITE_BEGIN(ModifierTestComposite)
TEST_COMPOSITE_ADD_UNIT(ModifierTest);
TEST_COMPOSITE_ADD_UNIT(ModifierManagerTest);
TEST_COMPOSITE_END

#endif
