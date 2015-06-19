#ifndef DDMD_STRESS_AUTO_CORRELATION_H
#define DDMD_STRESS_AUTO_CORRELATION_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/analyzers/AutoCorrAnalyzer.h>
#include <util/space/Tensor.h>

namespace DdMd
{

   using namespace Util;

   /**
   * Compute the shear stress autocorrelation function.
   *
   * This class uses a hierarchichal block averaging algorithm 
   * to calculate the shear stress autocorrelation function. 
   * The final *.dat file  function contains values of the 
   * quantity \f$ k_{B}T G(t) \f$, where \f$G(t)\f$ is the linear
   * shear stress modulus, i.e., the response of shear stress to 
   * a hypothetical infinitesimal step shear strain, per unit 
   * applied step strain. This quantity is calculated from the
   * autocorrelation of the traceless symmetric part of the total
   * stress, using symmetry relations appropriate for an isotropic
   * fluid.
   *
   * \ingroup DdMd_Analyzer_Module
   */
   class StressAutoCorrelation : public AutoCorrAnalyzer<Tensor, double>
   {
   
   public:
   
      /**
      * Constructor.
      *
      * \param simulation parent Simulation object. 
      */
      StressAutoCorrelation(Simulation& simulation);
   
      /**
      * Destructor.
      */
      virtual ~StressAutoCorrelation()
      {} 

      using AutoCorrAnalyzer<Tensor, double>::readParameters;
      using AutoCorrAnalyzer<Tensor, double>::loadParameters;
      using AutoCorrAnalyzer<Tensor, double>::save;
      using AutoCorrAnalyzer<Tensor, double>::clear;
      using AutoCorrAnalyzer<Tensor, double>::setup;
      using AutoCorrAnalyzer<Tensor, double>::sample;
      using AutoCorrAnalyzer<Tensor, double>::output;

   protected:

      virtual void computeData();
      virtual Tensor data();

   };

}
#endif 
