/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McAverageAnalyzer.h"
#include <mcMd/analyzers/base/AverageAnalyzer.tpp>
#include <mcMd/mcSimulation/McSystem.h>

namespace McMd
{

   /*
   * Constructor.
   */
   McAverageAnalyzer::McAverageAnalyzer(McSystem& system) 
    : AverageAnalyzer<McSystem>(system)
   {  setClassName("McAverageAnalyzer"); }

}
