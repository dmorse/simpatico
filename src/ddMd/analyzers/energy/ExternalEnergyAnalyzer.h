#ifndef DDMD_EXTERNAL_ENERGY_ANALYZER_H
#define DDMD_EXTERNAL_ENERGY_ANALYZER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/analyzers/AverageAnalyzer.h>

namespace DdMd
{

   using namespace Util;

   /**
   * Output and evaluate average of external energy.
   *
   * \sa \ref ddMd_analyzer_ExternalEnergyAnalyzer_page "param file format"
   *
   * \ingroup DdMd_Analyzer_Energy_Module
   */
   class ExternalEnergyAnalyzer : public AverageAnalyzer
   {
   
   public:
   
      /**
      * Constructor.
      *
      * \param simulation parent Simulation object. 
      */
      ExternalEnergyAnalyzer(Simulation& simulation);
   
      /**
      * Destructor.
      */
      virtual ~ExternalEnergyAnalyzer(); 
   
   protected:

      /**
      * Function to compute value.
      *
      * Call on all processors.
      */
      virtual void compute();

      /**
      * Current value, set by compute function.
      *
      * Call only on master.
      */
      virtual double value();
   
   };

}
#endif 
