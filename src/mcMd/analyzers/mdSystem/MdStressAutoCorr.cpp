/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/


#include "MdStressAutoCorr.h"        // class header

namespace McMd
{

   /* 
   * Constructor.
   */
   MdStressAutoCorr::MdStressAutoCorr(MdSystem& system)
    : StressAutoCorr<MdSystem>(system)
   {  setClassName("MdStressAutoCorr"); }

   /* 
   * Destructor.
   */
   MdStressAutoCorr::~MdStressAutoCorr()
   {}

   /* 
   * Compute stress.
   */
   void MdStressAutoCorr::computeStress(Tensor& stress)
   {
      Tensor virial;
      Tensor kinetic;
      system().computeVirialStress(virial);
      system().computeKineticStress(kinetic);
      stress.add(virial, kinetic);
   }

}
