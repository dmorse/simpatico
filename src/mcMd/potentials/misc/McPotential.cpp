/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/potentials/misc/McPotential.h>

namespace McMd
{

   /**
   * Constructor.
   */
   McPotential::McPotential(bool createsStress)
    : ParamComposite(),
      EnergyCalculator(),
      StressCalculator(createsStress)
   {  setClassName("McPotential"); }

   /**
   * Destructor (does nothing)
   */
   McPotential::~McPotential()
   {}

} 
