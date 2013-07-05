#ifndef MCMD_COULOMB_POTENTIAL_CPP
#define MCMD_COULOMB_POTENTIAL_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <stdlib.h>
//#include <util/space/Vector.h>
//#include <util/space/Tensor.h>
//#include <util/math/Constants.h>
#include <mcMd/simulation/System.h>
#include <mcMd/simulation/Simulation.h>
#include "CoulombPotential.h" 
#include "EwaldCoulombPair.h" 

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   CoulombPotential::CoulombPotential(System& system)
    : simulationPtr_(&system.simulation()),
      systemPtr_(&system),
      boundaryPtr_(&system.boundary()),
      atomTypesPtr_(&system.simulation().atomTypes()),
      pairInteractionPtr_(0),
      epsilon_(1.0),
      alpha_(1.0),
      rCutoff_(1.0),
      isInitialized_(false)
   {  setClassName("CoulombPotential"); }

   /*
   * Destructor (does nothing)
   */
   CoulombPotential::~CoulombPotential()
   {}

   /*
   * Create association with a EwaldCoulombPair.
   */
   void CoulombPotential::setPairInteraction(EwaldCoulombPair& pairInteraction)
   {
      pairInteractionPtr_ = &pairInteraction;
      pairInteractionPtr_->set(epsilon_, alpha_, rCutoff_);
   }

   /*
   * Read parameters and initialize.
   */
   void CoulombPotential::readParameters(std::istream& in)
   {
      read<double>(in, "epsilon", epsilon_);
      read<double>(in, "alpha",   alpha_);
      read<double>(in, "rCutoff", rCutoff_);
      if (pairInteractionPtr_) {
         pairInteractionPtr_->set(epsilon_, alpha_, rCutoff_);
      }
      isInitialized_ = true;
   }
 
} 
#endif
