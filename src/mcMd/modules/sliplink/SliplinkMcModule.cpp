#ifndef SLIPLINK_MC_MODULE_CPP
#define SLIPLINK_MC_MODULE_CPP

#include "SliplinkMcModule.h"
#include <mcMd/mcSimulation/McSimulation.h>
#include <mcMd/mcSimulation/McSystem.h>

// Custom factory classes
#include "analyzers/SliplinkMcAnalyzerFactory.h"
#include "mcMoves/SliplinkMcMoveFactory.h"

namespace McMd 
{

   using namespace Util;
   
   SliplinkMcModule::SliplinkMcModule(McSimulation& sim)
    : McModule(sim)
   {
      analyzerFactoryPtr_ 
          = new SliplinkMcAnalyzerFactory(simulation(), system());
      mcMoveFactoryPtr_     
          = new SliplinkMcMoveFactory(simulation(), system());
   }

   SliplinkMcModule::~SliplinkMcModule()
   {
      delete analyzerFactoryPtr_;
      delete mcMoveFactoryPtr_;
   }

   void SliplinkMcModule::addFactories()
   {
      simulation().analyzerFactory().addSubfactory(*analyzerFactoryPtr_);
      simulation().mcMoveFactory().addSubfactory(*mcMoveFactoryPtr_);
   }

}
#endif
