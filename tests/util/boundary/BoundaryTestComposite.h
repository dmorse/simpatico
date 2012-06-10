#ifndef BOUNDARY_TEST_COMPOSITE_H
#define BOUNDARY_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "OrthorhombicBoundaryTest.h"

TEST_COMPOSITE_BEGIN(BoundaryTestComposite)
TEST_COMPOSITE_ADD_UNIT(OrthorhombicBoundaryTest);
TEST_COMPOSITE_END

#endif
