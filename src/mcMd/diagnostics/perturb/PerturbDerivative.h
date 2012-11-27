#ifdef MCMD_PERTURB
#ifndef MCMD_PERTURB_DERIVATIVE_H
#define MCMD_PERTURB_DERIVATIVE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
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
   * PerturbDerivative returns average value of Perturbation::derivative().
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
      virtual void readParameters(std::istream& in);

      /**
      * Load internal state from an archive.
      *
      * \param ar input/loading archive
      */
      virtual void loadParameters(Serializable::IArchive &ar);
   
      /**
      * Save internal state to an archive.
      *
      * \param ar output/saving archive
      */
      virtual void save(Serializable::OArchive &ar);

      /**
      * Serialize to/from an archive. 
      *
      * \param ar      saving or loading archive
      * \param version archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

      /** 
      * Clear accumulator.
      */
      virtual void setup();

      /* 
      * Evaluate derivative and add to ensemble. 
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
      
      /// Index of perturbation parameter associated with derivative.
      int parameterIndex_;
      
      /// Has readParam been called?
      bool    isInitialized_;
      
   };

   /*
   * Serialize to/from an archive. 
   */
   template <class Archive>
   void PerturbDerivative::serialize(Archive& ar, const unsigned int version)
   {
      Diagnostic::serialize(ar, version);
      ar & nSamplePerBlock_;
      ar & parameterIndex_;
      ar & accumulator_;
   }

}
#endif   // ifndef PERTURB_DERIVATIVE_H
#endif   // ifdef  MCMD_PERTURB
