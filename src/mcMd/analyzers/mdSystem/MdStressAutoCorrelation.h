#ifndef MCMD_MD_STRESS_AUTOCORRELATION_H
#define MCMD_MD_STRESS_AUTOCORRELATION_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/system/StressAutoCorrelation.h>  // base class template
#include <mcMd/mdSimulation/MdSystem.h>                   // base template parameter

namespace McMd
{

   /**
   * Analyzer to calculate average isotropic pressure.
   */
   class MdStressAutoCorrelation : public StressAutoCorrelation<MdSystem>
   {
   public:

      /**
      * Constructor.
      *
      * \param system parent McSystem
      */
      MdStressAutoCorrelation(MdSystem& system);

      /**
      * Destructor.
      */
      ~MdStressAutoCorrelation();
   };

}
#endif 
