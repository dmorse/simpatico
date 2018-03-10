/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/potentials/special/SpecialPotential.h>

namespace McMd
{

   /**
   * Constructor.
   */
   SpecialPotential::SpecialPotential(bool createsStress)
    : ParamComposite(),
      EnergyCalculator(),
      StressCalculator(createsStress)
   {  setClassName("SpecialPotential"); }

   /**
   * Destructor (does nothing)
   */
   SpecialPotential::~SpecialPotential()
   {}

} 
