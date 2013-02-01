#ifndef SYSTEM_DIAGNOSTIC_FACTORY
#define SYSTEM_DIAGNOSTIC_FACTORY

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "SystemDiagnosticFactory.h" // Class header

// Diagnostics for any System (Mc or Md)
#include <mcMd/diagnostics/simulation/LogProgress.h>
#include "DumpConfig.h"
#include "RDF.h"
#include "StructureFactorP.h"
#include "StructureFactorPGrid.h"
#include "StructureFactor.h"
#include "StructureFactorGrid.h"
#include "IntraStructureFactor.h"
#include "VanHove.h"
#ifdef UTIL_MPI
#include "MigratingVanHove.h"
#endif
#include "RadiusGyration.h"
#include "BlockRadiusGyration.h"
#include "BondLengthDist.h"
#include "CompositionProfile.h"
#include "AtomMSD.h"
#include "ComMSD.h"
#include "IntraPairAutoCorr.h"
#include "RingRouseAutoCorr.h"
#include "LinearRouseAutoCorr.h"
#include "BoundaryAverage.h"

#ifdef MCMD_PERTURB
#ifdef UTIL_MPI
#include <mcMd/diagnostics/perturb/BennettsMethod.h>
#endif
#include <mcMd/diagnostics/perturb/PerturbDerivative.h>
#endif

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   SystemDiagnosticFactory::SystemDiagnosticFactory(Simulation& simulation, 
                                            System& system)
    : simulationPtr_(&simulation),
      systemPtr_(&system)
   {}

   /* 
   * Return a pointer to a instance of Species subclass className.
   */
   Diagnostic* SystemDiagnosticFactory::factory(const std::string &className) const
   {
      Diagnostic* ptr = 0;

      if (className == "LogProgress") {
         ptr = new LogProgress();
      } else
      if (className == "DumpConfig") {
         ptr = new DumpConfig(system());
      } else
      if (className == "RDF") {
         ptr = new RDF(system());
      } else 
      if (className == "StructureFactorP") {
         ptr = new StructureFactorP(system());
      } else 
      if (className == "StructureFactorPGrid") {
         ptr = new StructureFactorPGrid(system());
      } else 
      if (className == "StructureFactor") {
         ptr = new StructureFactor(system());
      } else 
      if (className == "StructureFactorGrid") {
         ptr = new StructureFactorGrid(system());
      } else 
      if (className == "IntraStructureFactor") {
         ptr = new IntraStructureFactor(system());
      } else 
      if (className == "VanHove") {
         ptr = new VanHove(system());
      } else 
      if (className == "BoundaryAverage") {
         ptr = new BoundaryAverage(system());
      } else 
      if (className == "RadiusGyration") {
         ptr = new RadiusGyration(system());
      } else 
      if (className == "BlockRadiusGyration") {
         ptr = new BlockRadiusGyration(system());
      } else 
      if (className == "BondLengthDist") {
         ptr = new BondLengthDist(system());
      } else 
      if (className == "CompositionProfile") {
         ptr = new CompositionProfile(system());
      } else 
      if (className == "AtomMSD") {
         ptr = new AtomMSD(system());
      } else
      if (className == "ComMSD") {
         ptr = new ComMSD(system());
      } else
      if (className == "IntraPairAutoCorr") {
         ptr = new IntraPairAutoCorr(system());
      } else
      if (className == "LinearRouseAutoCorr") {
         ptr = new LinearRouseAutoCorr(system());
      } else
      if (className == "RingRouseAutoCorr") {
         ptr = new RingRouseAutoCorr(system());
      } 

      #ifdef MCMD_PERTURB
      else 
      if (className == "PerturbDerivative") {
         ptr = new PerturbDerivative(system());
      }
      #ifdef UTIL_MPI
      else 
      if (className == "BennettsMethod") {
         ptr = new BennettsMethod(system());
      } else 
      if (className == "MigratingVanHove") {
         ptr = new MigratingVanHove(system());
      }

      #endif
      #endif

      return ptr;
   }

}

#endif
