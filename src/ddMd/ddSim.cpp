#ifndef DDMD_DDSIM_CPP
#define DDMD_DDSIM_CPP

#include <ddMd/system/System.h>

int main()
{
   DdMd::System system;
   system.readParam(std::cin); 
   system.readCommands();
}
#endif
