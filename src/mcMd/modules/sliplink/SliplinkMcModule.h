#ifndef SLIPLINK_MC_MODULE_H
#define SLIPLINK_MC_MODULE_H

#include <mcMd/mcSimulation/McModule.h>

namespace Util
{  template <class T> class Factory; }

namespace McMd 
{

   using namespace Util;

   class McMove;
   class Analyzer;

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
      * Add Analyzer and McMove sub-factories.
      */
      virtual void addFactories();
   
   private:
   
      Factory<Analyzer>* analyzerFactoryPtr_;

      Factory<McMove>*     mcMoveFactoryPtr_;
   
   };

}
#endif
