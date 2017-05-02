#ifndef HOOMD_MC_MODULE_CPP
#define HOOMD_MC_MODULE_CPP

#include "HoomdMcModule.h"
#include <mcMd/mcSimulation/McSimulation.h>
#include <mcMd/mcSimulation/McSystem.h>

// Custom factory classes
#include "analyzers/HoomdAnalyzerFactory.h"
#include "mcMoves/HoomdMcMoveFactory.h"
#include "potentials/pair/HoomdPairFactory.h"

#ifdef MCMD_PERTURB
#include "perturbation/HoomdMcPerturbationFactory.h"
#endif

#ifdef SIMP_EXTERNAL
#include "potentials/external/HoomdExternalFactory.h"
#include <mcMd/potentials/external/ExternalFactory.h>
#endif 

namespace McMd 
{

   using namespace Util;
   
   HoomdMcModule::HoomdMcModule(McSimulation& simulation)
    : McModule(simulation)
   {
      analyzerFactoryPtr_
          = new HoomdAnalyzerFactory(simulation, system());
      pairFactoryPtr_ 
          = new HoomdPairFactory();
      mcMoveFactoryPtr_     
          = new HoomdMcMoveFactory(simulation, system());
      #ifdef MCMD_PERTURB
      perturbationFactoryPtr_
          = new HoomdMcPerturbationFactory(system());
      #endif
      #ifdef SIMP_EXTERNAL
      externalFactoryPtr_ 
          = new HoomdExternalFactory(system());
      #endif
   }

   HoomdMcModule::~HoomdMcModule()
   {
      delete analyzerFactoryPtr_;
      delete pairFactoryPtr_;
      delete mcMoveFactoryPtr_;
      #ifdef MCMD_PERTURB
      delete perturbationFactoryPtr_;
      #endif
   }

   void HoomdMcModule::addFactories()
   {
      system().simulation().analyzerFactory().addSubfactory(*analyzerFactoryPtr_);
      #ifndef SIMP_NOPAIR
      system().pairFactory().addSubfactory(*pairFactoryPtr_);
      #endif
      simulation().mcMoveFactory().addSubfactory(*mcMoveFactoryPtr_);
      #ifdef MCMD_PERTURB
      if (system().expectPerturbation()) {
         system().perturbationFactory().addSubfactory(
            *perturbationFactoryPtr_);
      }
      #endif
      #ifdef SIMP_EXTERNAL
      system().externalFactory().addSubfactory(*externalFactoryPtr_);
      #endif
   }

}
#endif
