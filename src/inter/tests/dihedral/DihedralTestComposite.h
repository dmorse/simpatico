#ifndef DIHEDRAL_TEST_COMPOSITE_H
#define DIHEDRAL_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "TorsionForceTest.h"
#include "CosineDihedralTest.h"
#include "MultiHarmonicDihedralTest.h"

TEST_COMPOSITE_BEGIN(DihedralTestComposite)
TEST_COMPOSITE_ADD_UNIT(TorsionForceTest);
TEST_COMPOSITE_ADD_UNIT(CosineDihedralTest);
TEST_COMPOSITE_ADD_UNIT(MultiHarmonicDihedralTest);
TEST_COMPOSITE_END

#endif
