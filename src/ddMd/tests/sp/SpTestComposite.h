#ifndef DDMD_SP_TEST_COMPOSITE_CPP
#define DDMD_SP_TEST_COMPOSITE_CPP


#include "chemistry/SpSpeciesTest.h"
#include "storage/SpAtomStorageTest.h"
#include "storage/SpConfigurationTest.h"
#include "processor/ProcessorTest.h"
#include "configIos/DdMdSpConfigIoTest.h"

#include <test/CompositeTestRunner.h>

TEST_COMPOSITE_BEGIN(SpTestComposite)
addChild(new TEST_RUNNER(SpSpeciesTest), "chemistry/");
addChild(new TEST_RUNNER(SpSpeciesTest), "chemistry/");
addChild(new TEST_RUNNER(SpAtomStorageTest), "storage/");
addChild(new TEST_RUNNER(SpConfigurationTest), "storage/");
addChild(new TEST_RUNNER(ProcessorTest), "processor/");
//addChild(new TEST_RUNNER(DdMdSpConfigIoTest), "configIos/");
TEST_COMPOSITE_END

#endif
