#ifndef DDMD_TEST_CPP
#define DDMD_TEST_CPP

/*
* This program runs all unit tests in the DdMd directory.
*/ 

#include "storage/AtomStorageTest.h"
#include "storage/BondStorageTest.h"
#include "configIos/ConfigIoTest.h"
#include "communicate/CommunicateTestComposite.h"
#include "neighbor/NeighborTestComposite.h"
#include "simulation/SimulationTest.h"

#ifdef UTIL_MPI
#ifndef TEST_MPI
#define TEST_MPI
#endif
#endif

#include <test/CompositeTestRunner.h>

TEST_COMPOSITE_BEGIN(DdMdNsTestComposite)
addChild(new TEST_RUNNER(AtomStorageTest), "storage/");
addChild(new TEST_RUNNER(BondStorageTest), "storage/");
addChild(new TEST_RUNNER(ConfigIoTest), "configIos/");
addChild(new CommunicateTestComposite, "communicate/");
addChild(new NeighborTestComposite, "neighbor/");
addChild(new TEST_RUNNER(SimulationTest), "simulation/");
TEST_COMPOSITE_END

int main() {

   #ifdef UTIL_MPI 
   MPI::Init();
   Util::IntVector::commitMpiType();
   Util::Vector::commitMpiType();

   DdMdNsTestComposite runner;
   runner.run();

   MPI::Finalize();
   #endif

}
#endif
