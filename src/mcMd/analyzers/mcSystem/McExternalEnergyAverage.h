#ifdef SIMP_EXTERNAL
#ifndef MCMD_MC_EXTERNAL_ENERGY_AVERAGE_H
#define MCMD_MC_EXTERNAL_ENERGY_AVERAGE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/util/AverageAnalyzer.h> // base class template
#include <mcMd/mcSimulation/McSystem.h>          // base template parameter

#include <cstdio> 

namespace McMd
{

   using namespace Util;

   /**
   * McExternalEnergyAverage averages of total external energy.
   *
   * See \ref mcMd_analyzer_McExternalEnergyAverage_page "here" for 
   * the parameter file format and any other user documentation.
   *
   * \ingroup McMd_Analyzer_Mc_Module
   */
   class McExternalEnergyAverage : public AverageAnalyzer<McSystem>
   {
   
   public:

      /**   
      * Constructor.
      *
      * \param system parent McSystem
      */
      McExternalEnergyAverage(McSystem& system);

      /* 
      * Evaluate external energy per particle, and add to ensemble. 
      *
      * \param iStep step counter index
      */
      virtual void sample(long iStep);
   
      using AverageAnalyzer<McSystem>::readParameters;
      using AverageAnalyzer<McSystem>::loadParameters;
      using AverageAnalyzer<McSystem>::save;
      using AverageAnalyzer<McSystem>::setup;
      using AverageAnalyzer<McSystem>::output;

   protected:
 
      using AverageAnalyzer<McSystem>::outputFile_;
      using AverageAnalyzer<McSystem>::nSamplePerBlock_;
      using AverageAnalyzer<McSystem>::accumulator_;

   };

}
#endif 
#endif 
