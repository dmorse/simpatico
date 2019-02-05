#include "ModifierTestComposite.h"

using namespace Util;
using namespace DdMd;

int main() 
{

   #ifdef UTIL_MPI 
   #ifdef TEST_MPI
   MPI_Init(&argc, &argv);
   DdMd::Modifier::initStatic();
   #endif
   #endif

   ModifierTestComposite runner;
   runner.run();

   #ifdef UTIL_MPI 
   #ifdef TEST_MPI
   MPI_Finalize();
   #endif
   #endif

}
