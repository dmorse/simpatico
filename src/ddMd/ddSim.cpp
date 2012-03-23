#ifndef DDMD_DDSIM_CPP
#define DDMD_DDSIM_CPP

#include <ddMd/system/System.h>

int main()
{
   #ifdef UTIL_MPI
   MPI::Init();
   DdMd::System system(MPI::COMM_WORLD);
   #endif
   system.readParam(std::cin); 
   system.readCommands();
   MPI::Finalize();
}
#endif
