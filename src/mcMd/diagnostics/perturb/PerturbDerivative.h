#ifdef MCMD_PERTURB
#ifndef MCMD_PERTURB_DERIVATIVE_H
#define MCMD_PERTURB_DERIVATIVE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/diagnostics/SystemDiagnostic.h>  // base class template
#include <mcMd/simulation/System.h>             // base template parameter
#include <util/accumulators/Average.h>          // member

#include <cstdio> 

namespace McMd
{

   using namespace Util;

   /**
   * PerturbDerivative average return value of Perturbation::derivative().
   *
   * \ingroup McMd_Diagnostic_Module
   */
   class PerturbDerivative : public SystemDiagnostic<System>
   {
   
   public:

      /**   
      * Constructor.
      */
      PerturbDerivative(System& system);

      /**
      * Read parameters and initialize.
      */
      virtual void readParam(std::istream& in);

      /* 
      * Evaluate energy per particle, and add to ensemble. 
      */
      virtual void sample(long iStep);
   
      /**
      * Output results at end of simulation.
      */
      virtual void output();

   private:

      /// Output file stream
      std::ofstream outputFile_;

      /// Average object - statistical accumulator
      Average  accumulator_;

      /// Number of samples per block average output.
      int nSamplePerBlock_;
      
      int parameterIndex_;

   };

}
#endif   // ifndef PERTURB_DERIVATIVE_H
#endif   // ifdef  MCMD_PERTURB
