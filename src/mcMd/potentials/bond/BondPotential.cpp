/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/potentials/bond/BondPotential.h>

namespace McMd
{

   /**
   * Constructor.
   */
   BondPotential::BondPotential()
    : ParamComposite()
   {  setClassName("BondPotential"); }

   /**
   * Destructor (does nothing)
   */
   BondPotential::~BondPotential()
   {}

} 
