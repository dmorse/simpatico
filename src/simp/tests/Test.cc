/*
* This program runs all unit tests in the simp directory.
*/ 

#ifdef  UTIL_MPI
#ifndef TEST_MPI
#define TEST_MPI
#endif
#endif

#include "interaction/InteractionTestComposite.h"
#include "species/SpeciesTestComposite.h"
#include "boundary/BoundaryTestComposite.h"
#include "analysis/AnalysisTestComposite.h"
#include "crystal/CrystalTestComposite.h"
#include <test/CompositeTestRunner.h>

using namespace Simp;

TEST_COMPOSITE_BEGIN(SimpNsTestComposite)
addChild(new InteractionTestComposite, "interaction/");
addChild(new SpeciesTestComposite, "species/");
addChild(new BoundaryTestComposite, "boundary/");
addChild(new AnalysisTestComposite, "analysis/");
addChild(new CrystalTestComposite, "crystal/");
TEST_COMPOSITE_END


int main(int argc, char* argv[])
{
   try {

      #ifdef UTIL_MPI
      MPI::Init();
      Vector::commitMpiType();
      IntVector::commitMpiType();
      #endif 
   
      SimpNsTestComposite runner;
   
      if (argc > 2) {
         UTIL_THROW("Too many arguments");
      }
      if (argc == 2) {
         runner.addFilePrefix(argv[1]);
      }

      // Run all unit tests
      int failures = runner.run();
   
      #ifdef UTIL_MPI
      MPI::Finalize();
      #endif 

      return (failures != 0);

   } catch (...) {

      std::cerr << "Uncaught exception in src/simp/tests/Test.cc" << std::endl;
      return 1;
   }
}
