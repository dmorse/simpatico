#ifndef ANGLE_TEST_COMPOSITE_H
#define ANGLE_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "BendForceTest.h"
#include "CosineAngleTest.h"
#include "CosineSqAngleTest.h"
#include "HarmonicAngleTest.h"

TEST_COMPOSITE_BEGIN(AngleTestComposite)
TEST_COMPOSITE_ADD_UNIT(BendForceTest);
TEST_COMPOSITE_ADD_UNIT(CosineAngleTest);
TEST_COMPOSITE_ADD_UNIT(CosineSqAngleTest);
TEST_COMPOSITE_ADD_UNIT(HarmonicAngleTest);
TEST_COMPOSITE_END

#endif
