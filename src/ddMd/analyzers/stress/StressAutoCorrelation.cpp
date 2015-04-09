/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/analyzers/stress/StressAutoCorrelation.h>
#include <ddMd/analyzers/AutoCorrAnalyzer.tpp>

#include <sstream>

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   StressAutoCorrelation::StressAutoCorrelation(Simulation& simulation) 
    : AutoCorrAnalyzer<Tensor, double>(simulation)
   {  setClassName("StressAutoCorrelation"); }


   /*
   * Sample the stress tensor.
   */
   void StressAutoCorrelation::computeData() 
   {  
      Simulation& sim = simulation();
      sim.computeVirialStress();
      sim.computeKineticStress();
   }

   /*
   * Sample the stress tensor.
   */
   Tensor StressAutoCorrelation::data() 
   {  
      Tensor virial  = simulation().virialStress();
      Tensor kinetic = simulation().kineticStress();
      Tensor stress;
      stress.add(virial, kinetic);

      // Remove trace
      double pressure = 0.0;
      int i, j;
      for (i = 0; i < Dimension; ++i) {
         pressure += stress(i,i);
      }
      pressure = pressure/double(Dimension);
      for (i = 0; i < Dimension; ++i) {
         stress(i,i) -= pressure;
      }
   
      double factor = sqrt(simulation().boundary().volume()/10.0);
      for (i = 0; i < Dimension; ++i) {
         for (j = 0; j < Dimension; ++j) {
            stress(i,j) *= factor;
         }
      }

      return stress;
   }

}
