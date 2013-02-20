#ifndef MCMD_MC_STRESS_AUTO_CORR_H
#define MCMD_MC_STRESS_AUTO_CORR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/diagnostics/SystemDiagnostic.h>   // base class template
#include <mcMd/mcSimulation/McSystem.h>          // base class template parameter
#include <util/accumulators/AutoCorr.h>          // member template
#include <util/space/Tensor.h>                   // template parameter

namespace McMd
{

   using namespace Util;


   /**
   * Shear stress autocorrelation function.
   * \ingroup McMd_Diagnostic_Module
   */
   class McStressAutoCorr : public SystemDiagnostic<McSystem>
   {
   
   public:
  
      /**
      * Constructor.
      *
      * \param system reference to parent System
      */
      McStressAutoCorr(McSystem &system);
  
      /** 
      * Read parameters from file.
      *
      * Input format:
      *
      *   - int    interval        sampling interval
      *   - string outputFileName  output file name
      *   - int    capacity        capacity of array of previous values.
      *
      * \param in input parameter stream
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
      * Set number of molecules and clear accumulator.
      */
      virtual void setup();
   
      /** 
      * Evaluate separation vector for all chains, add to ensemble.
      *
      * \param iStep counter for number of steps
      */
      virtual void sample(long iStep);
   
      /**
      * Output results after simulation is completed.
      */
      virtual void output();

   private:
   
      /// Output file stream
      std::ofstream outputFile_;

      /// Statistical accumulator
      AutoCorr<Tensor, double> accumulator_;

      /// Maximum length of sequence in AutoCorr.
      int capacity_;
   
      // Has readParam been called?
      bool isInitialized_;

   };

   /*
   * Serialize to/from an archive. 
   */
   template <class Archive>
   void McStressAutoCorr::serialize(Archive& ar, const unsigned int version)
   {
      Diagnostic::serialize(ar, version);
      ar & capacity_;
      ar & accumulator_;
   }

}
#endif
