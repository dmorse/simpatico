#ifndef UTIL_TEST_CPP
#define UTIL_TEST_CPP

/*
* This program runs all unit tests in the util directory.
*/ 

#include "containers/ContainersTestComposite.h"
#include "archives/ArchiveTestComposite.h"
#include "param/serial/ParamTestComposite.h"
#include "space/SpaceTestComposite.h"
#include "crystal/CrystalTestComposite.h"
#include "format/FormatTest.h"
#include "random/RandomTest.h"

#include <test/CompositeTestRunner.h>

TEST_COMPOSITE_BEGIN(UtilNsTestComposite)
addChild(new ContainersTestComposite, "containers/");
addChild(new ParamTestComposite, "param/serial/");
addChild(new ArchiveTestComposite, "archives/");
addChild(new SpaceTestComposite, "space/");
addChild(new CrystalTestComposite, "crystal/");
addChild(new TEST_RUNNER(FormatTest), "format/");
addChild(new TEST_RUNNER(RandomTest), "random/");
TEST_COMPOSITE_END

int main() {

   UtilNsTestComposite runner;
   runner.run();

}
#endif
