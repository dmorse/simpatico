/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/potentials/misc/MdPotential.h>

namespace McMd
{

   /**
   * Constructor.
   */
   MdPotential::MdPotential(bool hasStress)
    : ParamComposite(),
      hasStress_(hasStress)
   {  setClassName("MdPotential"); }

   /**
   * Destructor (does nothing)
   */
   MdPotential::~MdPotential()
   {}

} 
