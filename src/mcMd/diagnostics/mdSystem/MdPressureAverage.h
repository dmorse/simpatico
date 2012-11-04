#ifndef MCMD_MD_PRESSURE_AVERAGE_H
#define MCMD_MD_PRESSURE_AVERAGE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/diagnostics/system/PressureAverage.h>  // base class template
#include <mcMd/mdSimulation/MdSystem.h>               // base template parameter

namespace McMd
{

   /**
   * Calculate average scalar pressure for an MdSystem.
   */
   class MdPressureAverage : public PressureAverage<MdSystem>
   {
   public:

      /*
      * Constructor
      *
      * \param system parent MdSystem
      */
      MdPressureAverage(MdSystem& system);

   };

}
#endif 
