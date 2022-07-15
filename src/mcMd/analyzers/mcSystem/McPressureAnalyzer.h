#ifndef MCMD_MC_PRESSURE_ANALYZER_H
#define MCMD_MC_PRESSURE_ANALYZER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/base/AverageAnalyzer.h> // base class template
#include <mcMd/mcSimulation/McSystem.h>          // base class parameter

namespace McMd
{

   using namespace Util;

   /**
   * McPressureAnalyzer evaluates average pressure in MC simulation.
   *
   * \ingroup McMd_Analyzer_McMd_Module
   */
   class McPressureAnalyzer : public AverageAnalyzer<McSystem>
   {

   public:

      /**   
      * Constructor.
      */
      McPressureAnalyzer(McSystem& system);

      /**
      * Compute pressure.
      */
      virtual void compute();

   };

}
#endif
