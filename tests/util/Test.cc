#ifndef UTIL_TEST_CPP
#define UTIL_TEST_CPP

/*
* This program runs all unit tests in the util directory.
*/ 

#include "accumulators/unit/AccumulatorTestComposite.h"
#include "archives/ArchiveTestComposite.h"
#include "boundary/BoundaryTestComposite.h"
#include "containers/ContainersTestComposite.h"
#include "crystal/CrystalTestComposite.h"
#include "format/FormatTest.h"
#include "param/serial/ParamTestComposite.h"
#include "random/RandomTest.h"
#include "space/SpaceTestComposite.h"

#include <test/CompositeTestRunner.h>

TEST_COMPOSITE_BEGIN(UtilNsTestComposite)
addChild(new AccumulatorTestComposite, "accumulators/unit/");
addChild(new ArchiveTestComposite, "archives/");
addChild(new BoundaryTestComposite, "boundary/");
addChild(new ContainersTestComposite, "containers/");
addChild(new CrystalTestComposite, "crystal/");
addChild(new TEST_RUNNER(FormatTest), "format/");
addChild(new ParamTestComposite, "param/serial/");
addChild(new TEST_RUNNER(RandomTest), "random/");
addChild(new SpaceTestComposite, "space/");
TEST_COMPOSITE_END

int main() {

   UtilNsTestComposite runner;
   runner.run();

}
#endif
