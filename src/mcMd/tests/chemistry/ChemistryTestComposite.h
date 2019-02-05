#ifndef MCMD_CHEMISTRY_TEST_COMPOSITE_H
#define MCMD_CHEMISTRY_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "AtomTypeTest.h"
#include "MaskPolicyTest.h"
#ifdef UTIL_MPI
#include "MpiChemistryTest.h"
#endif

TEST_COMPOSITE_BEGIN(ChemistryTestComposite)
TEST_COMPOSITE_ADD_UNIT(AtomTypeTest);
TEST_COMPOSITE_ADD_UNIT(MaskPolicyTest);
#ifdef UTIL_MPI
TEST_COMPOSITE_ADD_UNIT(MpiChemistryTest);
#endif
TEST_COMPOSITE_END

#endif
