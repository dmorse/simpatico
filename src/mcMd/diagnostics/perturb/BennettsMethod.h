#ifdef MCMD_PERTURB
#ifdef UTIL_MPI
#ifndef MCMD_BENNETTS_METHOD_H
#define MCMD_BENNETTS_METHOD_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/diagnostics/SystemDiagnostic.h>  // base class template
#include <mcMd/simulation/System.h>             // base template parameter
#include <util/accumulators/Average.h>          // member
#ifdef UTIL_MPI
#include <util/containers/DArray.h>
#include <util/mpi/MpiSendRecv.h>
#include <util/mpi/MpiLogger.h>
#endif
#include <mcMd/simulation/Simulation.h>
#include <util/space/Vector.h>          // Util namespace
#include <util/global.h>
#include <util/param/ParamComposite.h>

#include <fstream>
#include <cstdio> 

namespace McMd
{

   using namespace Util;

   /**
   * Bennet's Method gives an estimate of the free energy difference.
   *
   * \ingroup McMd_Diagnostic_Module
   */
   class BennettsMethod : public SystemDiagnostic<System>
   {
   
   public:

      /**   
      * Constructor.
      */
      BennettsMethod(System& system);

      /**
      * Read parameters and initialize.
      */
      virtual void readParam(std::istream& in);

      /* 
      * Evaluate energy per particle, and add to ensemble. 
      */
      virtual void sample(long iStep);
      virtual void initialize();
   
      /**
      * Output results at end of simulation.
      */
      virtual void output();
   
   protected:
      /// Value of the shift constant for the associated System.
      double shift_;
      double lowerShift_;

      #if UTIL_MPI

      // Values of shift constants for all processors.
      DArray<double> shifts_;

      #endif
     
   private:
      /// Tempering variable.
      DArray<double> myParam_;
      DArray<double> lowerParam_;
      DArray<double> upperParam_;
      double myDerivative_;

      /// System reference.
      System* systemPtr_;

      /// Get the communicator in the simulation.
      MPI::Intracomm* communicatorPtr_;
      double myArg_;
      double lowerArg_;
      double myFermi_;
      double upperFermi_;
      double lowerFermi_;

      /// Number of processors.
      int   nProcs_;

      /// Current processor's rank.
      int   myId_;

      /// Active neighboring (partner) replica's rank.
      int   lowerId_;
      int   upperId_;

      /// Number of simulation steps between subsequent actions.
      long interval_;

      /// Tags for exchanging parameters.
      static const int TagDerivative[2];
      static const int TagFermi[2];

      /// Output file stream
      std::ofstream outputFile_;

      /// Average object - statistical accumulator
      Average  myAccumulator_;
      Average  upperAccumulator_;

      /// Number of samples per block average output.
      int nSamplePerBlock_;

      int nParameters_;
  
      virtual void analyze();
   };

}
#endif   // ifndef BENNETTS_METHOD_H
#endif   // ifdef  UTIL_MPI
#endif   // ifdef  UTIL_PERTURB
