#ifndef DDMD_TEST_CPP
#define DDMD_TEST_CPP

/*
* This program runs all unit tests in the DdMd directory.
*/ 

#ifdef UTIL_MPI
#ifndef TEST_MPI
#define TEST_MPI
#endif
#endif

#include "chemistry/ChemistryTestComposite.h"
#include "storage/StorageTestComposite.h"
#include "configIos/ConfigIoTest.h"
#include "communicate/CommunicateTestComposite.h"
#include "neighbor/NeighborTestComposite.h"
#include "simulation/SimulationTest.h"
#ifdef DDMD_MODIFIERS
#include "modifiers/ModifierTestComposite.h"
#endif

#include <test/CompositeTestRunner.h>

TEST_COMPOSITE_BEGIN(DdMdNsTestComposite)
addChild(new ChemistryTestComposite, "chemistry/");
addChild(new StorageTestComposite, "storage/");
addChild(new NeighborTestComposite, "neighbor/");
#ifdef DDMD_MODIFIERS
addChild(new ModifierTestComposite, "modifiers/");
#endif
#ifdef TEST_MPI
addChild(new TEST_RUNNER(ConfigIoTest), "configIos/");
addChild(new CommunicateTestComposite, "communicate/");
addChild(new TEST_RUNNER(SimulationTest), "simulation/");
#endif
TEST_COMPOSITE_END

int main(int argc, char* argv[])
{
   // Precondition
   if (argc > 2) {
      UTIL_THROW("Too many arguments");
   }

   #ifdef TEST_MPI 
   MPI::Init();
   Util::IntVector::commitMpiType();
   Util::Vector::commitMpiType();
   #endif

   DdMdNsTestComposite runner;
   if (argc == 2) {
      runner.addFilePrefix(argv[1]);
   }
   runner.run();

   #ifdef TEST_MPI 
   MPI::Finalize();
   #endif

}
#endif
