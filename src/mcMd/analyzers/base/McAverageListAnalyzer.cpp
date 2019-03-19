/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McAverageListAnalyzer.h"
#include <mcMd/analyzers/base/AverageListAnalyzer.tpp>
#include <mcMd/mcSimulation/McSystem.h>

namespace McMd
{

   /*
   * Constructor.
   */
   McAverageListAnalyzer::McAverageListAnalyzer(McSystem& system) 
    : AverageListAnalyzer<McSystem>(system)
   {  setClassName("McAverageListAnalyzer"); }

}
