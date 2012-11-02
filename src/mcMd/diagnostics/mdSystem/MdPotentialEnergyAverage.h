#ifndef MCMD_MD_POTENTIAL_ENERGY_AVERAGE_H
#define MCMD_MD_POTENTIAL_ENERGY_AVERAGE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/diagnostics/SystemDiagnostic.h> // base class template
#include <mcMd/mdSimulation/MdSystem.h>        // class template parameter
#include <util/accumulators/Average.h>         // member

namespace McMd
{

   using namespace Util;

   /**
   * MdPotentialEnergyAverage averages of total potential energy.
   *
   * \ingroup McMd_Diagnostic_Module
   */
   class MdPotentialEnergyAverage : public SystemDiagnostic<MdSystem>
   {
   
   public:

      /**   
      * Constructor.
      */
      MdPotentialEnergyAverage(MdSystem& system);

      /**
      * Read parameters and initialize.
      */
      virtual void readParameters(std::istream& in);

      /**
      * Load state from an archive.
      *
      * \param ar loading (input) archive.
      */
      virtual void loadParameters(Serializable::IArchive& ar);

      /**
      * Save state to an archive.
      *
      * \param ar saving (output) archive.
      */
      virtual void save(Serializable::OArchive& ar);

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

      // Has readParam been called?
      bool isInitialized_;

   };

   /*
   * Serialize to/from an archive. 
   */
   template <class Archive>
   void MdPotentialEnergyAverage::serialize(Archive& ar, const unsigned int version)
   {
      Diagnostic::serialize(ar, version);
      ar & nSamplePerBlock_;
      ar & accumulator_;
   }

}
#endif 
