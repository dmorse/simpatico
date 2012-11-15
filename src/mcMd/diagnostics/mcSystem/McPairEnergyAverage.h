#ifndef INTER_NOPAIR
#ifndef MCMD_MC_PAIR_ENERGY_AVERAGE_H
#define MCMD_MC_PAIR_ENERGY_AVERAGE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/diagnostics/SystemDiagnostic.h>   // base class template
#include <mcMd/mcSimulation/McSystem.h>          // base template parameter
#include <util/accumulators/Average.h>           // member
#include <mcMd/neighbor/CellList.h>              // member
#include <mcMd/diagnostics/util/PairSelector.h>  // member

#include <cstdio> 

namespace McMd
{

   using namespace Util;

   /**
   * McPairEnergyAverage averages of total potential energy.
   *
   * \ingroup McMd_Diagnostic_Module
   */
   class McPairEnergyAverage : public SystemDiagnostic<McSystem>
   {
   
   public:

      /**   
      * Constructor.
      *
      * \param system parent McSystem
      */
      McPairEnergyAverage(McSystem& system);

      /**
      * Read parameters and initialize.
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream& in);

      /**
      * Load state from an archive.
      *
      * \param ar loading (input) archive
      */
      virtual void loadParameters(Serializable::IArchive& ar);

      /**
      * Save state to an archive.
      *
      * \param ar saving (output) archive
      */
      virtual void save(Serializable::OArchive& ar);
  
      /**
      * Serialize to/from an archive. 
      * 
      * \param ar      archive
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

      /// List of neighbors 
      mutable CellList::NeighborArray neighbors_;

      /// Selector to determine which pairs to include.
      PairSelector  selector_;

      // Has readParam been called?
      bool isInitialized_;

   };

   /*
   * Serialize to/from an archive. 
   */
   template <class Archive>
   void McPairEnergyAverage::serialize(Archive& ar, const unsigned int version)
   {
      Diagnostic::serialize(ar, version);
      ar & nSamplePerBlock_;
      ar & selector_;
      ar & accumulator_;
   }

}
#endif 
#endif 
