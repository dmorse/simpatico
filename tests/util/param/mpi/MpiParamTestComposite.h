#ifndef MPI_PARAM_TEST_COMPOSITE_H
#define MPI_PARAM_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "MpiParamTest.h"
#include "MpiParamCompositeTest.h"
#include "MpiManagerTest.h"

TEST_COMPOSITE_BEGIN(MpiParamTestComposite)
TEST_COMPOSITE_ADD_UNIT(MpiParamTest);
TEST_COMPOSITE_ADD_UNIT(MpiParamCompositeTest);
TEST_COMPOSITE_ADD_UNIT(MpiManagerTest);
TEST_COMPOSITE_END

#endif
