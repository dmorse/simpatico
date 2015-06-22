#ifndef DDMD_STRESS_ANALYZER_H
#define DDMD_STRESS_ANALYZER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/analyzers/SymmTensorAverageAnalyzer.h>

namespace DdMd
{

   using namespace Util;

   /**
   * Output and evaluate average of stress tensor.
   *
   * \sa \ref ddMd_analyzer_StressAnalyzer_page "param file format"
   *
   * \ingroup DdMd_Analyzer_Module
   */
   class StressAnalyzer : public SymmTensorAverageAnalyzer
   {

   public:

      /**
      * Constructor.
      *
      * \param simulation parent Simulation object.
      */
      StressAnalyzer(Simulation& simulation);

      /**
      * Destructor.
      */
      virtual ~StressAnalyzer();

   protected:

      /**
      * Compute current value of stress tensor.
      *
      * Call on all processors.
      */
      virtual void compute();

      /**
      * Get current value, set by compute function.
      *
      * Call only on master.
      */
      virtual Tensor value();

   };

}
#endif
