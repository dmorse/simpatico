#ifndef UTIL_TEST_CPP
#define UTIL_TEST_CPP

/*
* This program runs all unit tests in the util directory.
*/ 

#include <test/CompositeTestRunner.h>

#include <util/containers/tests/ContainersTestComposite.h>
#include <util/archives/tests/ArchiveTestComposite.h>
#include <util/param/tests/serial/ParamTestComposite.h>
#include <util/space/tests/SpaceTestComposite.h>
#include <util/crystal/tests/CrystalTestComposite.h>
#include <util/format/tests/FormatTest.h>
#include <util/random/tests/RandomTest.h>

TEST_COMPOSITE_BEGIN(UtilNsTestComposite)
addChild(new ContainersTestComposite, "containers/tests/");
addChild(new ParamTestComposite, "param/tests/serial/");
addChild(new ArchiveTestComposite, "archives/tests/");
addChild(new SpaceTestComposite, "space/tests/");
addChild(new CrystalTestComposite, "crystal/tests/");
addChild(new TEST_RUNNER(FormatTest), "format/tests/");
addChild(new TEST_RUNNER(RandomTest), "random/tests/");
TEST_COMPOSITE_END

int main() {

   UtilNsTestComposite runner;
   runner.run();

}
#endif
