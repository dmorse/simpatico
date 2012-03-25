#ifndef DDMD_EXTERNAL_POTENTIAL_CPP
#define DDMD_EXTERNAL_POTENTIAL_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "ExternalPotential.h"
#include <ddMd/simulation/Simulation.h>
#include <ddMd/storage/AtomStorage.h>
#include <ddMd/storage/AtomIterator.h>
#include <ddMd/communicate/Domain.h>
#include <util/space/Vector.h>
#include <util/global.h>

namespace DdMd
{
   using namespace Util;

   /*
   * Constructor.
   */
   ExternalPotential::ExternalPotential(Simulation& simulation)
    : simulationPtr_(&simulation),
      boundaryPtr_(&simulation.boundary()),
      domainPtr_(&simulation.domain()),
      storagePtr_(&simulation.atomStorage())
   {}

   /*
   * Constructor (for unit testing).
   */
   ExternalPotential::ExternalPotential(Boundary& boundary, Domain& domain,
                                AtomStorage& storage)
    : boundaryPtr_(&boundary),
      domainPtr_(&domain),
      storagePtr_(&storage)
   {} 

   /*
   * Destructor.
   */
   ExternalPotential::~ExternalPotential()
   {}

}
#endif
