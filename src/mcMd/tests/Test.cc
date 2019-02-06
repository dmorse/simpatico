/*
* This program runs all unit tests in the mcMd directory.
*/

#ifdef UTIL_MPI
#define TEST_MPI
#endif

#include <test/CompositeTestRunner.h>

#include "chemistry/ChemistryTestComposite.h"
#include "neighbor/NeighborTestComposite.h"
#include "simulation/SimulationTestComposite.h"
#include "mcSimulation/McSimulationTest.h"
#include "mdSimulation/MdSimulationTest.h"
#include "util/misc/initStatic.h"

// Define a class McMdNsTestComposite
TEST_COMPOSITE_BEGIN(McMdNsTestComposite)
addChild(new ChemistryTestComposite, "chemistry/");
addChild(new NeighborTestComposite, "neighbor/");
addChild(new SimulationTestComposite, "simulation/");
addChild(new TEST_RUNNER(McSimulationTest), "mcSimulation/");
addChild(new TEST_RUNNER(MdSimulationTest), "mdSimulation/");
#ifdef UTIL_MPI
addChild(new TEST_RUNNER(MpiChemistryTest), "chemistry/");
#endif
TEST_COMPOSITE_END

int main(int argc, char* argv[])
{
   Util::initStatic();
   #ifdef UTIL_MPI
   MPI_Init(&argc, &argv);
   McMd::commitMpiTypes();
   #endif

   try {
      if (argc > 2) {
         UTIL_THROW("Too many arguments");
      }
   
      McMdNsTestComposite runner;
      if (argc == 2) {
         runner.addFilePrefix(argv[1]);
       }
   
      // Run all unit test methods
      int failures = runner.run();

      #ifdef UTIL_MPI
      MPI_Finalize();
      #endif

      return (failures != 0);

   } catch (...) {

      //std::cerr << "Uncaught exception in mcMd/tests/Test.cc" << std::endl;
      UTIL_THROW("Uncaught exception in mcMd/tests/Test.cc");
      //return 1;

   }
}
