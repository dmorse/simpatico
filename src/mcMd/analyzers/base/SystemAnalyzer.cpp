/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SystemAnalyzer.h"
#include <mcMd/simulation/System.h>
#include <mcMd/mcSimulation/McSystem.h>
#include <mcMd/mdSimulation/MdSystem.h>

namespace McMd
{

   template class SystemAnalyzer<System>;
   template class SystemAnalyzer<McSystem>;
   template class SystemAnalyzer<MdSystem>;

}
