#ifndef SLIPLINK_MC_MODULE_CPP
#define SLIPLINK_MC_MODULE_CPP

#include "SliplinkMcModule.h"
#include <mcMd/mcSimulation/McSimulation.h>
#include <mcMd/mcSimulation/McSystem.h>

// Custom factory classes
#include "diagnostics/SliplinkMcDiagnosticFactory.h"
#include "mcMoves/SliplinkMcMoveFactory.h"

namespace McMd 
{

   using namespace Util;
   
   SliplinkMcModule::SliplinkMcModule(McSimulation& sim)
    : McModule(sim)
   {
      diagnosticFactoryPtr_ 
          = new SliplinkMcDiagnosticFactory(simulation(), system());
      mcMoveFactoryPtr_     
          = new SliplinkMcMoveFactory(simulation(), system());
   }

   SliplinkMcModule::~SliplinkMcModule()
   {
      delete diagnosticFactoryPtr_;
      delete mcMoveFactoryPtr_;
   }

   void SliplinkMcModule::addFactories()
   {
      simulation().diagnosticFactory().addSubfactory(*diagnosticFactoryPtr_);
      simulation().mcMoveFactory().addSubfactory(*mcMoveFactoryPtr_);
   }

}
#endif
