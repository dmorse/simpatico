#ifndef CHEMISTRY_TEST_COMPOSITE_H
#define CHEMISTRY_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "AtomTypeTest.h"
#include "SpeciesGroupTest.h"

TEST_COMPOSITE_BEGIN(ChemistryTestComposite)
TEST_COMPOSITE_ADD_UNIT(AtomTypeTest);
TEST_COMPOSITE_ADD_UNIT(SpeciesGroupTest);
TEST_COMPOSITE_END

#endif
