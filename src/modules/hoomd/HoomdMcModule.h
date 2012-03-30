#ifndef HOOMD_MC_MODULE_H
#define HOOMD_MC_MODULE_H

#include <mcMd/mcSimulation/McModule.h>

#include <util/param/Factory.h>

#ifdef MCMD_PERTURB
#include <mcMd/perturb/Perturbation.h>
#endif

namespace McMd 
{

   using namespace Util;

   class McMove;
   class PairFactory;

   #ifdef INTER_EXTERNAL
   class ExternalFactory;
   #endif

   /**
   * Module for slip link simulation. 
   */
   class HoomdMcModule : public McModule
   {

   public:

      /**
      * Constructor.  
      */
      HoomdMcModule(McSimulation& simulation);

      /**
      * Destructor.  
      */
      ~HoomdMcModule();

      /**
      * Add specialized Pair and McMove factories.
      */
      virtual void addFactories();
   
   private:
   
      Factory<McMove>*     mcMoveFactoryPtr_;

      PairFactory*         pairFactoryPtr_;

      #ifdef MCMD_PERTURB
      Factory<Perturbation>* perturbationFactoryPtr_;
      #endif   

      #ifdef INTER_EXTERNAL
      ExternalFactory*     externalFactoryPtr_;
      #endif
   };

}
#endif
