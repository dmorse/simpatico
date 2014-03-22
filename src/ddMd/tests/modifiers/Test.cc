#include "ModifierTestComposite.h"

using namespace Util;
using namespace DdMd;

int main() 
{

   #ifdef UTIL_MPI 
   #ifdef TEST_MPI
   MPI::Init();
   DdMd::Modifier::initStatic();
   #endif
   #endif

   ModifierTestComposite runner;
   runner.run();

   #ifdef UTIL_MPI 
   #ifdef TEST_MPI
   MPI::Finalize();
   #endif
   #endif

}
