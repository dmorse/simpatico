#ifndef SIMP_ANALYSIS_TEST_COMPOSITE_H
#define SIMP_ANALYSIS_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "AverageMixInTest.h"

TEST_COMPOSITE_BEGIN(AnalysisTestComposite)
TEST_COMPOSITE_ADD_UNIT(AverageMixInTest);
TEST_COMPOSITE_END

#endif
