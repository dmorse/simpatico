#ifdef MCMD_PERTURB
#ifdef UTIL_MPI
#ifndef BENNETTS_METHOD_H
#define BENNETTS_METHOD_H

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
   * \ingroup Diagnostic_Module
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
      
      /// Shift factor (the logarithm of the difference between free 
      /// energies of replicas n and n+1) associated with replica n.
      double shift_;
      
      /// Shift factor (the logarithm of the difference between free 
      /// energies of replicas n-1 and n) associated with replica n.
      double lowerShift_;

      #if UTIL_MPI

      // Values of shift constants for all processors.
      DArray<double> shifts_;

      #endif
     
   private:
      
      /// Get the communicator in the simulation.
      MPI::Intracomm* communicatorPtr_;
      
      /// Current processor's rank.
      int   myId_;
      
      /// Number of processors.
      int   nProcs_;

      /// Active neighboring (lower partner) replica's rank.
      int   lowerId_;
      
      /// Active neighboring (upper partner) replica's rank.
      int   upperId_;

      /// Number of perturbation parameters.
      int nParameters_;
     
      /// Average object - statistical accumulator
      Average  myAccumulator_;
      
      /// Average object - statistical accumulator
      Average  upperAccumulator_;

      /// Number of samples per block average output.
      int nSamplePerBlock_;
 
      /// Tempering parameter.
      DArray<double> myParam_;
      
      /// Lower partner's tempering parameter.
      DArray<double> lowerParam_;
      
      /// Upper partner's tempering parameter.
      DArray<double> upperParam_;
      
      /// Derivative of hamiltonian with respect to a perturbation parameter variable.
      double myDerivative_;
      
      /// Argument representing the difference between free energy estimated from
      /// derivative calculations and the input shift factor (for replicas n and n+1).
      double myArg_;
      
      /// Argument representing the difference between free energy estimated from
      /// derivative calculations and the input shift factor (for replicas n and n-1).
      double lowerArg_;

      /// Fermi function of the argument myArg_.
      double myFermi_;

      /// Fermi function of the argument lowerArg_.
      double lowerFermi_;

      /// Fermi function received from upper replica.
      double upperFermi_;
      
      /// Output file stream
      std::ofstream outputFile_;

      /// Tags for exchanging parameters.
      static const int TagDerivative[2];
      static const int TagFermi[2];
  
      virtual void analyze();
   };

}
#endif   // ifndef BENNETTS_METHOD_H
#endif   // ifdef  UTIL_MPI
#endif   // ifdef  UTIL_PERTURB
