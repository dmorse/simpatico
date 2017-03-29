#include "SliplinkMcModule.h"

namespace McMd 
{

   using namespace Util;
   
   SliplinkMcModule::SliplinkMcModule(McSimulation& sim)
    : analyzerFactory_(sim, sim.system()),
      mcMoveFactory_(sim, sim.system())
   {
      sim.analyzerFactory().addSubfactory(analyzerFactory_);
      sim.mcMoveFactory().addSubfactory(mcMoveFactory_);
   }

}
