#ifndef MCMD_SLIPLINK_MC_MODULE_H
#define MCMD_SLIPLINK_MC_MODULE_H

// Custom factory classes
#include <mcMd/mcSimulation/McSimulation.h>
#include "analyzers/SliplinkMcAnalyzerFactory.h"
#include "mcMoves/SliplinkMcMoveFactory.h"

namespace McMd 
{

   using namespace Util;

   /**
   * Module for slip link simulation. 
   */
   class SliplinkMcModule 
   {

   public:

      /**
      * Constructor.  
      *
      * \param sim parent McSimulation.
      */
      SliplinkMcModule(McSimulation& sim);

   private:
   
      SliplinkMcAnalyzerFactory analyzerFactory_;
      SliplinkMcMoveFactory mcMoveFactory_;
   
   };

}
#endif
