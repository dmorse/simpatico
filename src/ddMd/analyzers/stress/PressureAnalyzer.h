#ifndef DDMD_PRESSURE_ANALYZER_H
#define DDMD_PRESSURE_ANALYZER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/analyzers/base/AverageAnalyzer.h>

namespace DdMd
{

   using namespace Util;

   /**
   * Sample and evaluate average of total pressure.
   *
   * \sa \ref ddMd_analyzer_PressureAnalyzer_page "param file format"
   *
   * \ingroup DdMd_Analyzer_Stress_Module
   */
   class PressureAnalyzer : public AverageAnalyzer
   {

   public:

      /**
      * Constructor.
      *
      * \param simulation  parent Simulation object.
      */
      PressureAnalyzer(Simulation& simulation);

      /**
      * Destructor.
      */
      virtual ~PressureAnalyzer();

   protected:

      /**
      * Compute current value of pressure.
      *
      * Call on all processors.
      */
      virtual void compute();

      /**
      * Get current value of pressure.
      *
      * Call only on master.
      */
      virtual double value();

   };

}
#endif
