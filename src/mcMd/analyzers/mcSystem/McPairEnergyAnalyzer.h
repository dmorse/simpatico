#ifndef SIMP_NOPAIR
#ifndef MCMD_MC_PAIR_ENERGY_ANALYZER_H
#define MCMD_MC_PAIR_ENERGY_ANALYZER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/base/AverageAnalyzer.h> // base class template
#include <mcMd/mcSimulation/McSystem.h>          // base class parameter
#include <util/accumulators/Average.h>           // member
#include <mcMd/neighbor/CellList.h>              // member
#include <mcMd/analyzers/util/PairSelector.h>    // member

namespace McMd
{

   using namespace Util;

   /**
   * McPairEnergyAnalyzer averages pair energy for one type of pair.
   *
   * See \ref mcMd_analyzer_McPairEnergyAnalyzer_page "here" for 
   * the parameter file format and any other user documentation.
   *
   * \ingroup McMd_Analyzer_Mc_Module
   */
   class McPairEnergyAnalyzer : public AverageAnalyzer<McSystem>
   {
   
   public:

      /**   
      * Constructor.
      *
      * \param system parent McSystem
      */
      McPairEnergyAnalyzer(McSystem& system);

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
  
      /* 
      * Evaluate energy per particle, and add to ensemble. 
      */
      virtual void compute();
   
   private:

      /// List of neighbors 
      mutable CellList::NeighborArray neighbors_;

      /// Selector to determine which pairs to include.
      PairSelector  selector_;

   };

}
#endif 
#endif 
