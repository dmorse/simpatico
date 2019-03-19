#ifndef MCMD_MC_AVERAGE_ANALYZER_H
#define MCMD_MC_AVERAGE_ANALYZER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/base/AverageListAnalyzer.h>  // base class template

namespace McMd
{

   class McSystem;

   /**
   * AverageListAnalyzer for MC simulations.
   *
   * \ingroup McMd_Analyzer_Base_Module
   */
   class McAverageListAnalyzer : public AverageListAnalyzer<McSystem>
   {
   
   public:

      /**   
      * Constructor.
      */
      McAverageListAnalyzer(McSystem& system);

   };

}
#endif 
