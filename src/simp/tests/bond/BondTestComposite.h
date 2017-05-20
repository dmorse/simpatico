#ifndef BOND_TEST_COMPOSITE_H
#define BOND_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "HarmonicBondTest.h"
#include "HarmonicL0BondTest.h"

TEST_COMPOSITE_BEGIN(BondTestComposite)
TEST_COMPOSITE_ADD_UNIT(HarmonicBondTest);
TEST_COMPOSITE_ADD_UNIT(HarmonicL0BondTest);
TEST_COMPOSITE_END

#endif
