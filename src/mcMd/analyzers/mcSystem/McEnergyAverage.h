#ifndef MCMD_MC_ENERGY_AVERAGE_H
#define MCMD_MC_ENERGY_AVERAGE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/base/McAverageAnalyzer.h>  // base class 

#include <cstdio> 

namespace McMd
{

   using namespace Util;

   /**
   * McEnergyAverage averages of total potential energy.
   *
   * See \ref mcMd_analyzer_McEnergyAverage_page "here" for the
   * parameter file format and any other user documentation.
   *
   * \ingroup McMd_Analyzer_Mc_Module
   */
   class McEnergyAverage : public McAverageAnalyzer
   {
   
   public:

      /**   
      * Constructor.
      */
      McEnergyAverage(McSystem& system);

   protected:
 
      virtual void compute();

   };

}
#endif 
