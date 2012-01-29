#ifndef SLIPLINK_MC_MODULE_H
#define SLIPLINK_MC_MODULE_H

#include <mcMd/mcSimulation/McModule.h>

namespace Util
{  template <class T> class Factory; }

namespace McMd 
{

   using namespace Util;

   class McMove;
   class Diagnostic;

   /**
   * Module for slip link simulation. 
   */
   class SliplinkMcModule : public McModule
   {

   public:

      /**
      * Constructor.  
      *
      * \param sim parent McSimulation.
      */
      SliplinkMcModule(McSimulation& sim);

      /**
      * Destructor.  
      */
      virtual ~SliplinkMcModule();

      /**
      * Add Diagnostic and McMove sub-factories.
      */
      virtual void addFactories();
   
   private:
   
      Factory<Diagnostic>* diagnosticFactoryPtr_;

      Factory<McMove>*     mcMoveFactoryPtr_;
   
   };

}
#endif
