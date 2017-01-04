#ifndef MCMD_MD_STRESS_AUTOCORR_H
#define MCMD_MD_STRESS_AUTOCORR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/system/StressAutoCorr.h> // base class template
#include <mcMd/mdSimulation/MdSystem.h>           // base template parameter

namespace McMd
{

   /**
   * Analyzer to calculate average isotropic pressure.
   *
   * See \ref mcMd_analyzer_MdStressAutoCorr_page "here" for the
   * parameter file format and any other user documentation.
   *
   * This analyzer computes the stress autocorrelation function
   * \f[
   *    C(t) = k_{B}T G(t)
   * \f]
   * for an isotropic liquid, where \f$G(t)\f$ is the shear 
   * stress relaxation modulus. This function is given by a sum
   * \f[
   *    C(t) = \frac{1}{10} \sum_{ij=0}^{2} 
   *           \langle V
   *           [\sigma_{ij}(t) - P(t)\delta_{ij}]
   *           [\sigma_{ij}(0) - P(0)\delta_{ij}]
   *           \rangle
   * \f]
   * where \f$V = \sqrt{V(t)V(0)}\f$ is the system volume, 
   * \f$\sigma_{ij}(t)\f$ is the total stress tensor (virial + kinetic), 
   * and 
   * \f[ 
   *    P(t) \equiv \sum_{i=0}^{2} \sigma_{ii}(t)/3
   * \f]
   * is the pressure.
   *
   * \ingroup McMd_Analyzer_Md_Module
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
      * \param stress computed stress tensor (on return).
      */ 
      void computeStress(Tensor& stress);

   };

}
#endif 
