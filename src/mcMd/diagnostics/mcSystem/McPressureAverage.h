#ifndef MC_PRESSURE_AVERAGE_H
#define MC_PRESSURE_AVERAGE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/diagnostics/system/PressureAverage.h>  // base class template
#include <mcMd/mcSimulation/McSystem.h>               // base template parameter

namespace McMd
{

   /**
   * Diagnostic to calculate average isotropic pressure.
   */
   class McPressureAverage : public PressureAverage<McSystem>
   {
   
   public:

      /// Constructor.
      McPressureAverage(McSystem& system);

   };

}
#endif 
