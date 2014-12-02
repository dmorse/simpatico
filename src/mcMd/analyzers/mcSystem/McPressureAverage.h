#ifndef MCMD_MC_PRESSURE_AVERAGE_H
#define MCMD_MC_PRESSURE_AVERAGE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/system/PressureAverage.h>  // base class template
#include <mcMd/mcSimulation/McSystem.h>               // base template parameter

namespace McMd
{

   /**
   * Analyzer to calculate average isotropic pressure.
   */
   class McPressureAverage : public PressureAverage<McSystem>
   {
   public:

      /**
      * Constructor.
      *
      * \param system parent McSystem
      */
      McPressureAverage(McSystem& system);

      /**
      * Destructor.
      */
      ~McPressureAverage();

      // Use inherited PressureAverage<McSystem>::loadParameters
      // Use inherited PressureAverage<McSystem>::save

   };

}
#endif 
