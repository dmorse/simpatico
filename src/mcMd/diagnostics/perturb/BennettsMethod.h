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
      *       
      * \param system reference to parent System object
      */
      BennettsMethod(System& system);

      /**
      * Read parameters from file.
      *  
      * \param in input parameter stream
      */
      virtual void readParam(std::istream& in);

      /** 
      * Clear accumulators.
      */
      virtual void initialize();

      /* 
      * Evaluate the Fermi functions and add to accumulators.
      *       
      * \param iStep step counter
      */
      virtual void sample(long iStep);
   
      /**
      * Output results at end of simulation.
      */
      virtual void output();
   
   protected:
  
      /// Value of the shift constant for the associated system.
      double shift_;
      
      /// Value of the shift constant for the lower replica system.
      double lowerShift_;

      #if UTIL_MPI

      // Values of shift constants for all replicas.
      DArray<double> shifts_;

      #endif
     
   private:

      /// Tags for exchanging parameters.
      static const int TagDerivative[2];
      
      /// Tags for exchanging parameters.
      static const int TagFermi[2];
      
      /// Get the communicator in the simulation.
      MPI::Intracomm* communicatorPtr_;
      
      /// Current processor's rank.
      int   myId_;

      /// Number of processors.
      int   nProcs_;
      
      /// Lower replica's rank.
      int   lowerId_;
      
      /// Upper replica's rank.
      int   upperId_;
      
      /// Number of perturbation parameters.
      int nParameters_;
      
      /// Tempering variable.
      DArray<double> myParam_;

      /// Tempering variable of lower replica.
      DArray<double> lowerParam_;

      /// Tempering variable of upper replica.
      DArray<double> upperParam_;

      /// Number of samples per block average output.
      int nSamplePerBlock_;

      /// Average object - statistical accumulator.
      Average  myAccumulator_;
      
      /// Average object associated with samples from upper replica.
      Average  upperAccumulator_;

      /// Used to store an argument associated with this system.
      double myArg_;
      
      /// Used to store an argument associated with system of lower replica.
      double lowerArg_;
      
      /// Fermi function of argument associated with this system.
      double myFermi_;
      
      /// Fermi function of argument associated with system of lower replica.
      double lowerFermi_;
      
      /// Fermi function of argument associated with system of upper replica.
      double upperFermi_;

      /// Output file stream.
      std::ofstream outputFile_;
      
      /// Number of simulation steps between subsequent actions.
      long interval_;
  
      /**
      * Average of the collected samples used to estimate free energy 
      * differences at end of simulation.
      */
      virtual void analyze();
   };

}
#endif   // ifndef BENNETTS_METHOD_H
#endif   // ifdef  UTIL_MPI
#endif   // ifdef  UTIL_PERTURB
