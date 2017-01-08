#ifndef MCMD_MD_PRESSURE_AVERAGE_H
#define MCMD_MD_PRESSURE_AVERAGE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/system/PressureAverage.h>  // base class template
#include <mcMd/mdSimulation/MdSystem.h>               // base template parameter

namespace McMd
{

   /**
   * Calculate average scalar pressure for an MdSystem.
   *
   * See \ref mcMd_analyzer_MdPressureAverage_page "here" for 
   * the parameter file format and any other user documentation.
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

      // using inherited PressureAverage<MdSystem>::loadParameters
      // using inherited PressureAverage<MdSystem>::save

   };

}
#endif 
