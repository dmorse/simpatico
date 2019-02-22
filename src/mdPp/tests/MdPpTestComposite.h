#ifndef MDPP_TEST_COMPOSITE_CPP
#define MDPP_TEST_COMPOSITE_CPP

#include "storage/AtomStorageTest.h"
//#include "storage/SpeciesStorageTest.h"
#include "storage/ConfigurationTest.h"
#include "processor/ProcessorTest.h"
#include "neighbor/CellTest.h"
#include "config/ConfigTestComposite.h"

#include <test/CompositeTestRunner.h>

TEST_COMPOSITE_BEGIN(MdPpTestComposite)
//addChild(new TEST_RUNNER(SpeciesTest), "chemistry/");
addChild(new TEST_RUNNER(AtomStorageTest), "storage/");
addChild(new TEST_RUNNER(ConfigurationTest), "storage/");
addChild(new TEST_RUNNER(ProcessorTest), "processor/");
addChild(new TEST_RUNNER(CellTest), "neighbor/");
addChild(new ConfigTestComposite, "config/");
TEST_COMPOSITE_END

#endif
