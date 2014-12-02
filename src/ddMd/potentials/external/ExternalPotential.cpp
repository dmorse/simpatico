/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ExternalPotential.h"
#include <ddMd/simulation/Simulation.h>
#include <ddMd/storage/AtomStorage.h>
#include <ddMd/storage/AtomIterator.h>
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
      storagePtr_(&simulation.atomStorage())
   {  setClassName("ExternalPotential"); }

   /*
   * Default constructor (for unit testing).
   */
   ExternalPotential::ExternalPotential()
    : boundaryPtr_(0),
      storagePtr_(0)
   {} 

   /*
   * Associate with related objects. (for unit testing).
   */
   void ExternalPotential::associate(Boundary& boundary, AtomStorage& storage)
   {
      boundaryPtr_ = &boundary;
      storagePtr_ = &storage;
   } 

   /*
   * Destructor.
   */
   ExternalPotential::~ExternalPotential()
   {}

}
