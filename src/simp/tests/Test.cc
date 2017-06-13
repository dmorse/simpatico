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
#include <test/CompositeTestRunner.h>

using namespace Simp;

TEST_COMPOSITE_BEGIN(SimpNsTestComposite)
addChild(new InteractionTestComposite, "interaction/");
addChild(new SpeciesTestComposite, "species/");
TEST_COMPOSITE_END


int main(int argc, char* argv[])
{
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
   runner.run();

   #ifdef UTIL_MPI
   MPI::Finalize();
   #endif 
}
