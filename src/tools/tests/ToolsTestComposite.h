#ifndef TOOLS_TEST_COMPOSITE_CPP
#define TOOLS_TEST_COMPOSITE_CPP

#include "chemistry/SpeciesTest.h"
#include "storage/AtomStorageTest.h"
#include "storage/ConfigurationTest.h"
#include "processor/ProcessorTest.h"
#include "neighbor/CellTest.h"
#include "config/DdMdConfigReaderTest.h"

#include <test/CompositeTestRunner.h>

TEST_COMPOSITE_BEGIN(ToolsTestComposite)
addChild(new TEST_RUNNER(SpeciesTest), "chemistry/");
addChild(new TEST_RUNNER(SpeciesTest), "chemistry/");
addChild(new TEST_RUNNER(AtomStorageTest), "storage/");
addChild(new TEST_RUNNER(ConfigurationTest), "storage/");
addChild(new TEST_RUNNER(ProcessorTest), "processor/");
addChild(new TEST_RUNNER(CellTest), "neighbor/");
addChild(new TEST_RUNNER(DdMdConfigReaderTest), "config/");
TEST_COMPOSITE_END

#endif
