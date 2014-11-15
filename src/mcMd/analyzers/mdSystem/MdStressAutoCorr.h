#ifndef MCMD_MD_STRESS_AUTOCORR_H
#define MCMD_MD_STRESS_AUTOCORR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/system/StressAutoCorr.h>  // base class template
#include <mcMd/mdSimulation/MdSystem.h>                   // base template parameter

namespace McMd
{

   /**
   * Analyzer to calculate average isotropic pressure.
   */
   class MdStressAutoCorr : public StressAutoCorr<MdSystem>
   {
   public:

      /**
      * Constructor.
      *
      * \param system parent McSystem
      */
      MdStressAutoCorr(MdSystem& system);

      /**
      * Destructor.
      */
      ~MdStressAutoCorr();

   protected:
  
      /**
      * Compute total stress tensor.
      * 
      * \param stress Stress tensor (on return).
      */ 
      void computeStress(Tensor& total);

   };

}
#endif 
