#ifndef MDPP_CONFIG_TEST_COMPOSITE_H
#define MDPP_CONFIG_TEST_COMPOSITE_H

#include "DdMdConfigReaderTest.h"
#include "SmpConfigReaderTest.h"
#include "HoomdConfigReaderTest.h"
#include "TypeMapTest.h"

#include <test/CompositeTestRunner.h>

TEST_COMPOSITE_BEGIN(ConfigTestComposite)
TEST_COMPOSITE_ADD_UNIT(DdMdConfigReaderTest);
TEST_COMPOSITE_ADD_UNIT(SmpConfigReaderTest);
TEST_COMPOSITE_ADD_UNIT(TypeMapTest);
TEST_COMPOSITE_ADD_UNIT(HoomdConfigReaderTest);
TEST_COMPOSITE_END

#endif
