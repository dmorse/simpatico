#ifndef PAIR_TEST_COMPOSITE_H
#define PAIR_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "LJPairTest.h"
#include "DpdPairTest.h"

TEST_COMPOSITE_BEGIN(PairTestComposite)
TEST_COMPOSITE_ADD_UNIT(LJPairTest);
TEST_COMPOSITE_ADD_UNIT(DpdPairTest);
TEST_COMPOSITE_END

#endif
