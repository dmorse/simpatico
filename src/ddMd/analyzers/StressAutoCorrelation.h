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
   * Calculate stress autocorrelation function.
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
