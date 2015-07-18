#ifndef MCMD_ANALYZERS_TEST_COMPOSITE_H
#define MCMD_ANALYZERS_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "ClusterTest.h"

TEST_COMPOSITE_BEGIN(AnalyzersTestComposite)
TEST_COMPOSITE_ADD_UNIT(ClusterTest);
TEST_COMPOSITE_END

#endif
