#ifndef DDMD_KINETIC_ENERGY_ANALYZER_H
#define DDMD_KINETIC_ENERGY_ANALYZER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/analyzers/AverageAnalyzer.h>
//#include <ddMd/simulation/Simulation.h>
//#include <util/accumulators/Average.h>

namespace DdMd
{

   using namespace Util;

   /**
   * Output and evaluate average of kinetic energy.
   *
   * \sa \ref ddMd_analyzer_KineticEnergyAnalyzer_page "param file format"
   *
   * \ingroup DdMd_Analyzer_Module
   */
   class KineticEnergyAnalyzer : public AverageAnalyzer
   {
   
   public:
   
      /**
      * Constructor.
      *
      * \param simulation parent Simulation object. 
      */
      KineticEnergyAnalyzer(Simulation& simulation);
   
      /**
      * Destructor.
      */
      virtual ~KineticEnergyAnalyzer(); 
   
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
