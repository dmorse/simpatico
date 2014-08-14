#ifndef SPAN_TEST_COMPOSITE_CPP
#define SPAN_TEST_COMPOSITE_CPP

#include "chemistry/SpeciesTest.h"
#include "storage/AtomStorageTest.h"
#include "storage/ConfigurationTest.h"
#include "processor/ProcessorTest.h"
#include "configIos/DdMdConfigIoTest.h"
#include "neighbor/CellTest.h"

#include <test/CompositeTestRunner.h>

TEST_COMPOSITE_BEGIN(SpAnTestComposite)
addChild(new TEST_RUNNER(SpeciesTest), "chemistry/");
addChild(new TEST_RUNNER(SpeciesTest), "chemistry/");
addChild(new TEST_RUNNER(AtomStorageTest), "storage/");
addChild(new TEST_RUNNER(ConfigurationTest), "storage/");
addChild(new TEST_RUNNER(ProcessorTest), "processor/");
addChild(new TEST_RUNNER(CellTest), "neighbor/");
addChild(new TEST_RUNNER(DdMdConfigIoTest), "configIos/");
TEST_COMPOSITE_END

#endif
