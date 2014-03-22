#ifndef DDMD_TESTS_CHEMISTRY_TEST_COMPOSITE_H
#define DDMD_TESTS_CHEMISTRY_TEST_COMPOSITE_H

#include "AtomTest.h"
#include "GroupTest.h"

#include <test/CompositeTestRunner.h>

TEST_COMPOSITE_BEGIN(ChemistryTestComposite)
TEST_COMPOSITE_ADD_UNIT(AtomTest);
TEST_COMPOSITE_ADD_UNIT(GroupTest);
TEST_COMPOSITE_END

#endif
