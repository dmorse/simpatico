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
   * Bennett's method estimates free energy difference between two states.
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
      virtual void readParameters(std::istream& in);

      /**
      * Load internal state from an archive.
      *
      * \param ar input/loading archive
      */
      virtual void loadParameters(Serializable::IArchive& ar);

      /**
      * Save internal state to an archive.
      *
      * \param ar output/saving archive
      */
      virtual void save(Serializable::OArchive& ar);

      /**
      * Serialize to/from an archive.  
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

      /** 
      * Clear accumulators.
      */
      virtual void setup();

      /* 
      * Evaluate Fermi functions and add to accumulators.
      *       
      * \param iStep step counter
      */
      virtual void sample(long iStep);
   
      /**
      * Output results at end of simulation.
      */
      virtual void output();
   
   protected:
  
      /// Value of shift constant for associated system.
      double shift_;
      
      /// Value of shift constant for lower replica system.
      double lowerShift_;

      #if UTIL_MPI

      // Values of shift constants for all replicas.
      DArray<double> shifts_;

      #endif
     
   private:

      /// Tags for exchanging derivatives.
      static const int TagDerivative[2];
      
      /// Tags for exchanging Fermi functions.
      static const int TagFermi[2];
      
      /// MPI communicator in the simulation.
      MPI::Intracomm* communicatorPtr_;
      
      /// Current processor's rank.
      int  myId_;

      /// Number of processors.
      int  nProcs_;
      
      /// Lower replica's rank.
      int  lowerId_;
      
      /// Upper replica's rank.
      int  upperId_;
      
      /// Number of perturbation parameters.
      int nParameter_;
      
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
      
      /// Average object associated with upper replica.
      Average  upperAccumulator_;

      /// Used to store an argument.
      double myArg_;
      
      /// Used to store an argument of lower replica.
      double lowerArg_;
      
      /// Fermi function of the argument variable.
      double myFermi_;
      
      /// Fermi function of argument of lower replica.
      double lowerFermi_;
      
      /// Fermi function of argument of upper replica.
      double upperFermi_;

      /// Output file stream.
      std::ofstream outputFile_;
      
      /// Sampling interval.
      long interval_;
 
      /// Has readParam been called?
      bool  isInitialized_;
      
      /**
      * Estimate free energy difference between two states (replicas) 
      * from accumulators of two states.
      */
      virtual void analyze();
   };

   /*
   * Serialize to/from an archive. 
   */
   template <class Archive>
   void BennettsMethod::serialize(Archive& ar, const unsigned int version)
   {
      Diagnostic::serialize(ar, version);
      ar & nSamplePerBlock_;
      #ifdef UTIL_MPI
      if (hasIoCommunicator()) {
         ar & shifts_;
      } else {
         ar & shift_;
      }
      #else
      ar & shift_;
      #endif
      ar & myAccumulator_; 
      ar & upperAccumulator_; 

      // Set in constructor
      //ar & myId_;
      //ar & nProcs_;
      //ar & lowerId_;
      //ar & upperId_;
      //ar & nParameter_;
      // ar & myParam_;
      
      ar & lowerParam_;
      ar & upperParam_;
      ar & myArg_;
      ar & lowerArg_;
      ar & myFermi_;
      ar & lowerFermi_;
      ar & upperFermi_;
   }


}
#endif   // ifndef BENNETTS_METHOD_H
#endif   // ifdef  UTIL_MPI
#endif   // ifdef  UTIL_PERTURB
